use nom::{Err::Incomplete, Needed};
use std::io::{BufRead, Read};

use crate::*;
/// Represents the state of the BAM Reader
///
/// Header => Next call to `read()` will parse BAM header
/// Reference => Next call to `read()` will parse references
/// Alignment => Next call to `read()` will parse an alignment record
/// Complete => Reader has been exhausted. Subsequent calls will only produce Complete.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BamReaderState {
    Header,
    Reference,
    Alignment,
    Complete,
}

/// A streaming BAM Reader
///
/// Accepts any source implementing BufRead
/// Assumes input is uncompressed so must be coupled with a blocked gzip reader for compressed data.
/// Can handle incomplete input for header and references, but must be able to read
/// an entire alignment block into memory at once.
pub struct BamReader<T>
where
    T: BufRead,
{
    inner: T,
    buffer: Vec<u8>,
    offset: usize,
    state: BamReaderState,
    pub header: Option<BamHeader>,
    pub references: Vec<BamReference>,
}

impl<T> BamReader<T>
where
    T: BufRead,
{
    pub fn new(handle: T) -> Self {
        BamReader {
            inner: handle,
            buffer: Vec::with_capacity(MAX_BLOCK_SIZE),
            offset: 0,
            state: BamReaderState::Header,
            header: None,
            references: Vec::with_capacity(1),
        }
    }

    fn get_slice(&self) -> &[u8] {
        &self.buffer[self.offset..]
    }

    fn read_header(&mut self) -> BamReaderState {
        self.read_to_buffer(8).unwrap();
        while self.header.is_none() {
            match parser::read_header(self.get_slice()) {
                Ok((_, res)) => {
                    self.header = Some(res);
                }
                Err(Incomplete(Needed::Size(s))) => {
                    self.read_to_buffer(u64::try_from(s.get()).unwrap())
                        .unwrap();
                }
                Err(e) => panic!("Unable to parse BAM header: {e}"),
            }
        }

        if self.header.as_ref().unwrap().n_ref > 0 {
            self.state = BamReaderState::Reference;
        } else {
            self.state = BamReaderState::Alignment;
        }
        self.buffer.clear();
        self.state
    }

    fn read_references(&mut self) -> BamReaderState {
        let n_ref = usize::try_from(self.header.as_ref().unwrap().n_ref).unwrap();
        self.references = Vec::with_capacity(n_ref);
        while self.references.len() < n_ref {
            match parser::read_reference(self.get_slice()) {
                Ok((i, bref)) => {
                    self.offset = self.buffer.len() - i.len();
                    self.references.push(bref);
                }
                Err(Incomplete(Needed::Size(s))) => {
                    self.read_to_buffer(u64::try_from(s.get()).unwrap())
                        .unwrap();
                }
                Err(e) => panic!("Malformed BAM reference: {e}"),
            }
        }
        self.buffer.clear();
        self.offset = 0;
        self.state = BamReaderState::Alignment;
        self.state
    }

    fn read_to_buffer(&mut self, amt: u64) -> Result<u64, std::io::Error> {
        std::io::copy(&mut self.inner.by_ref().take(amt), &mut self.buffer)
    }

    /// Attempt to read a full alignment block into buffer.
    ///
    /// Returns amount read if it is equal to the block_size field.
    /// If there is no more input to be read from inner reader, returns Ok(0), signaling EOF.
    fn read_block(&mut self) -> Result<u64, BamError> {
        match self.read_to_buffer(4u64) {
            Ok(4u64) => {}
            Ok(0) => return Ok(0),
            Ok(_) => return Err(BamError::EofError),
            Err(e) => return Err(BamError::IoError(e)),
        }
        match parser::block_size(self.get_slice()) {
            Ok((_, bsize)) => match self.read_to_buffer(u64::from(bsize)) {
                Ok(v) if v == u64::from(bsize) => Ok(v),
                Ok(_) => Err(BamError::EofError),
                Err(e) => Err(BamError::IoError(e)),
            },
            Err(_) => Err(BamError::ParseError),
        }
    }

    fn read_record(&mut self) -> Option<Result<Record, BamError>> {
        match self.state {
            BamReaderState::Alignment => {
                match self.read_block() {
                    Ok(0) => {
                        self.state = BamReaderState::Complete;
                        return None;
                    }
                    Err(e) => return Some(Err(e)),
                    _ => {}
                }
                match parser::read_alignment(self.get_slice(), &self.references) {
                    Ok((_, aln)) => {
                        self.buffer.clear();
                        Some(Ok(aln))
                    }
                    Err(_) => Some(Err(BamError::ParseError)),
                }
            }
            BamReaderState::Complete => None,
            BamReaderState::Header => {
                self.read_header();
                self.read_record()
            }
            BamReaderState::Reference => {
                self.read_references();
                self.read_record()
            }
        }
    }
}

impl<B> Iterator for BamReader<B>
where
    B: BufRead,
{
    type Item = Result<Record, BamError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_record()
    }
}
