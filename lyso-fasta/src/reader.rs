use crate::parser;
use crate::FastaError;
use crate::Record;
use nom::Err::Incomplete;
use std::io::BufRead;

const MAX_BUFFER_SIZE: usize = 10_000_000;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FastaReaderState {
    Reading,
    Complete,
    Failed,
}

pub struct FastaReader<T>
where
    T: BufRead,
{
    state: FastaReaderState,
    inner: T,
    buffer: Vec<u8>,
    offset: usize,
}

impl<T> FastaReader<T>
where
    T: BufRead,
{
    pub fn new(f: T) -> Self {
        let r = FastaReader {
            state: FastaReaderState::Reading,
            inner: f,
            buffer: Vec::with_capacity(MAX_BUFFER_SIZE),
            offset: 0,
        };
        r
    }

    /// Prevent internal buffer from growing infinitely.
    /// Does not shrink capacity under the assumption that
    /// reads in a fasta tend to be of similar length.
    #[inline]
    fn resize_buffer(&mut self) {
        {
            self.buffer.drain(0..self.offset);
        }
        self.offset = 0;
    }

    #[inline]
    fn get_slice(&self) -> &[u8] {
        &self.buffer[self.offset..]
    }

    #[inline]
    fn read_to_next_header(&mut self) -> Result<usize, std::io::Error> {
        self.inner.read_until(b'>', &mut self.buffer)
    }

    #[inline]
    pub fn read_record(&mut self) -> Option<Result<Record, FastaError>> {
        if self.state != FastaReaderState::Reading {
            return None;
        }
        match self.read_to_next_header() {
            Ok(0) if self.offset == self.buffer.len() => {
                self.state = FastaReaderState::Complete;
                return None;
            }
            Ok(_) => {}
            Err(e) => return Some(Err(FastaError::IoError(e))),
        }
        let mut res: Option<Result<Record, FastaError>> = None;
        while res.is_none() {
            match parser::parse_record(self.get_slice()) {
                Ok((i, (id, seq))) => {
                    self.offset = self.buffer.len() - i.len();
                    res = Some(Ok(Record { id, seq }))
                }
                Err(Incomplete(_)) => match self.read_to_next_header() {
                    Ok(0) => {
                        return Some(Err(FastaError::EofError));
                    }
                    Ok(_) => {}
                    Err(e) => return Some(Err(FastaError::IoError(e))),
                },
                Err(_) => {
                    self.state = FastaReaderState::Failed;
                    return Some(Err(FastaError::ParserError));
                }
            }
        }
        if self.offset > MAX_BUFFER_SIZE {
            self.resize_buffer();
        }
        res
    }
}

impl<T> Iterator for FastaReader<T>
where
    T: BufRead,
{
    type Item = Result<Record, FastaError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_record()
    }
}

#[cfg(test)]
mod tests {

    use crate::reader::FastaReader;
    use std::fs::File;
    use std::io::BufReader;

    const FA_PATH: &str = "../resources/test_data/test.fa";
    const BAD_FA_PATH: &str = "../resources/test_data/corrupt.fa";

    #[test]
    fn test_read_fa() {
        let f = File::open(FA_PATH).unwrap();
        let b = BufReader::new(f);
        let reader: FastaReader<BufReader<File>> = FastaReader::new(b);
        for r in reader {
            assert!(r.is_ok());
        }
    }

    #[test]
    #[should_panic]
    fn test_bad_fa_panics() {
        let f = File::open(BAD_FA_PATH).unwrap();
        let b = BufReader::new(f);
        let reader: FastaReader<BufReader<File>> = FastaReader::new(b);
        for r in reader {
            eprintln!("{:?}", r);
        }
    }

    #[test]
    fn test_get_fields() {
        let f = File::open(FA_PATH).unwrap();
        let b = BufReader::new(f);
        let mut reader: FastaReader<BufReader<File>> = FastaReader::new(b);
        let record = reader.next().unwrap();
        eprintln!("{}", record.as_ref().unwrap().seq);
        assert!(record.as_ref().unwrap().id == "SRR22092847.1.1");
        assert!(
            record.unwrap().seq
                == "GNTTAAAGCACATAAAGACAAATCGCTCCAGGGCAAAGNTTAAAGCACATAAAGACAAATCGCTCCAGGGCAAA"
        );
    }
}
