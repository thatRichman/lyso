use crate::parser;
use crate::FastaError;
use crate::Record;
use nom::Err::Incomplete;
use std::borrow::BorrowMut;
use std::io::BufRead;
use std::io::Read;

const MAX_BUFFER_SIZE: usize = 1_000_000;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FastaReaderState {
    Reading,
    Complete,
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
        FastaReader {
            state: FastaReaderState::Reading,
            inner: f,
            buffer: Vec::with_capacity(MAX_BUFFER_SIZE),
            offset: 0,
        }
    }

    fn maybe_resize_buffer(&mut self) {
        if self.buffer.len() > MAX_BUFFER_SIZE {
            eprintln!("Buffer is {}, resizing", self.buffer.len());
            self.buffer.drain(0..self.offset);
            self.buffer.shrink_to(MAX_BUFFER_SIZE);
            self.offset = 0;
        }
    }

    fn get_slice(&self) -> &[u8] {
        &self.buffer[self.offset..]
    }

    fn read_to_buffer(&mut self, n: u64) -> Result<usize, std::io::Error> {
        (&mut self.inner).take(n).read_to_end(&mut self.buffer)
    }

    fn read_to_next_header(&mut self) -> Result<usize, std::io::Error> {
        self.inner.read_until(b'>', &mut self.buffer)
    }

    pub fn read_record(&mut self) -> Option<Result<Record, FastaError>> {
        if self.state == FastaReaderState::Complete {
            return None;
        }
        match self.read_to_buffer(u64::try_from(MAX_BUFFER_SIZE).unwrap()) {
            Ok(0) if self.buffer.is_empty() => {
                self.state = FastaReaderState::Complete;
                return None;
            }
            Ok(v) => {
                eprintln!("read {} bytes", v)
            }
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
                    Ok(0) => return Some(Err(FastaError::EofError)),
                    Ok(v) => {
                        eprintln!("read {} bytes to next header", v)
                    }
                    Err(e) => return Some(Err(FastaError::IoError(e))),
                },
                Err(_) => return Some(Err(FastaError::ParserError)),
            }
        }
        self.maybe_resize_buffer();
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
    use std::path::PathBuf;

    fn init_path(s: &str) -> PathBuf {
        let mut test_data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        test_data_dir.push(s);
        test_data_dir
    }

    #[test]
    fn test_read_fa() {
        let fa_path = init_path("resources/test_data/test.fa");
        let f = File::open(fa_path).unwrap();
        let b = BufReader::new(f);
        let reader: FastaReader<BufReader<File>> = FastaReader::new(b);
        for r in reader {
            assert!(r.is_ok());
        }
    }

    #[test]
    #[should_panic]
    fn test_bad_fa_panics() {
        let fa_path = init_path("resources/test_data/corrupt.fa");
        let f = File::open(fa_path).unwrap();
        let b = BufReader::new(f);
        let reader: FastaReader<BufReader<File>> = FastaReader::new(b);
        for r in reader {
            eprintln!("{:?}", r);
        }
    }

    #[test]
    fn test_get_fields() {
        let fa_path = init_path("resources/test_data/test.fa");
        let f = File::open(fa_path).unwrap();
        let b = BufReader::new(f);
        let mut reader: FastaReader<BufReader<File>> = FastaReader::new(b);
        let record = reader.next().unwrap();

        assert!(record.unwrap().id == "SRR22092847.1.1");
        assert!(
            record.unwrap().seq
                == "GNTTAAAGCACATAAAGACAAATCGCTCCAGGGCAAAGNTTAAAGCACATAAAGACAAATCGCTCCAGGGCAAA"
        );
    }
}
