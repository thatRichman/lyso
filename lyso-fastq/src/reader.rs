use nom::Err::Incomplete;
use nom::Needed;
use std::io::BufRead;

use crate::parser;
use crate::{FastqError, Record};

const MAX_BUFFER_SIZE: usize = 10_000_000;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FastqReaderState {
    Reading,
    Complete,
    Failed,
}

#[derive(Debug)]
pub struct FastqReader<T> {
    state: FastqReaderState,
    inner: T,
    buffer: Vec<u8>,
    offset: usize,
}

impl<T> FastqReader<T>
where
    T: BufRead,
{
    pub fn new(f: T) -> Self {
        FastqReader {
            state: FastqReaderState::Reading,
            inner: f,
            buffer: Vec::with_capacity(MAX_BUFFER_SIZE),
            offset: 0,
        }
    }

    /// Prevent internal buffer from growing infinitely.
    /// Does not shrink capacity under the assumption that
    /// reads in a fastq tend to be of similar length.
    #[inline]
    fn resize_buffer(&mut self) {
        {
            self.buffer.drain(0..self.offset)
        };
        self.offset = 0;
    }

    #[inline]
    fn get_slice(&self) -> &[u8] {
        &self.buffer[self.offset..]
    }

    #[inline]
    /// FASTQ records are always 4 lines, so try to read that much
    fn read_to_buffer(&mut self) -> Result<usize, std::io::Error> {
        let mut amt = 0;
        for _ in 0..4 {
            amt += (&mut self.inner).read_until(b'\n', &mut self.buffer)?;
        }
        Ok(amt)
    }

    #[inline]
    pub fn read_record(&mut self) -> Option<Result<Record, FastqError>> {
        if self.state != FastqReaderState::Reading {
            return None;
        }
        match self.read_to_buffer() {
            Ok(0) if self.offset == self.buffer.len() => {
                self.state = FastqReaderState::Complete;
                return None;
            }
            Ok(_) => {}
            Err(e) => return Some(Err(FastqError::IoError(e))),
        }
        let mut res: Option<Result<Record, FastqError>> = None;
        while res.is_none() {
            match parser::parse_record(self.get_slice()) {
                Ok((i, (id, desc, seq, qual))) => {
                    res = Some(Ok(Record {
                        id: id.to_string(),
                        desc: desc.to_string(),
                        seq: seq.to_string(),
                        qual: qual.to_string(),
                    }));
                    self.offset = self.buffer.len() - i.len();
                }
                Err(Incomplete(Needed::Size(_))) => match self.read_to_buffer() {
                    Ok(0) => {
                        return Some(Err(FastqError::EofError));
                    }
                    Ok(_) => {}
                    Err(e) => return Some(Err(FastqError::IoError(e))),
                },
                Err(_) => {
                    return Some(Err(FastqError::ParseError));
                }
            }
        }
        if self.offset > MAX_BUFFER_SIZE {
            self.resize_buffer();
        }
        res
    }
}

impl<T> Iterator for FastqReader<T>
where
    T: BufRead,
{
    type Item = Result<Record, FastqError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_record()
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::PathBuf;

    fn init_path(s: &str) -> PathBuf {
        let mut test_data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        test_data_dir.push(s);
        test_data_dir
    }

    #[test]
    fn test_read_fq() {
        let fq_path = init_path("resources/test_data/test.fastq");
        let f = File::open(fq_path).unwrap();
        let b = BufReader::new(f);
        let reader = FastqReader::new(b);
        for record in reader {
            eprintln!("{:?}", record);
        }
    }

    #[test]
    #[should_panic]
    fn test_bad_fq_panics() {
        let fq_path = init_path("resources/test_data/corrupt.fastq");
        let f = File::open(fq_path).unwrap();
        let b = BufReader::new(f);
        let reader = FastqReader::new(b);
        for _ in reader {
            continue;
        }
    }

    #[test]
    fn test_get_fields() {
        let fq_path = init_path("resources/test_data/test.fastq");
        let f = File::open(fq_path).unwrap();
        let b = BufReader::new(f);
        let mut reader = FastqReader::new(b);
        let record = reader.next().expect("bad record!").expect("bad record!");

        assert!(record.id == "SRR22092847.1.1");
        assert!(record.desc == "1 length=37");
        assert!(record.qual == "F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
        assert!(record.seq == "GNTTAAAGCACATAAAGACAAATCGCTCCAGGGCAAA");
    }

    #[test]
    fn test_read_checked_fq() {
        let fq_path = init_path("resources/test_data/test.fastq");
        let f = File::open(fq_path).unwrap();
        let b = BufReader::new(f);
        let mut reader = FastqReader::new(b);
        assert!(reader.all(|x| x.is_ok()));
    }

    #[test]
    fn test_bad_fq_is_recoverable() {
        let fq_path = init_path("resources/test_data/trunc.fastq");
        let f = File::open(fq_path).unwrap();
        let b = BufReader::new(f);
        let reader = FastqReader::new(b);
        for record in reader {
            if let Some(r) = record.into() {
                match r {
                    Ok(_) => {}
                    Err(e) => eprintln!("{}", e),
                }
            }
        }
    }

    #[test]
    fn test_corrupt_fq_is_recoverable() {
        let fq_path = init_path("resources/test_data/corrupt.fastq");
        let f = File::open(fq_path).unwrap();
        let b = BufReader::new(f);
        let reader = FastqReader::new(b);
        for record in reader {
            if let Some(r) = record.into() {
                match r {
                    Ok(_) => {}
                    Err(e) => eprintln!("{}", e),
                }
            }
        }
    }
}
