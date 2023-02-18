use std::borrow::Cow;
use std::io::{BufRead, Seek};
use std::iter::Iterator;
use std::marker::PhantomData;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FastaError {
    #[error("Validation error")]
    ValidationError(&'static str),
    #[error("Unexpected end of file")]
    EofError,
    #[error("Missing id field")]
    MissingId,
    #[error("Missing sequence")]
    MissingSequenceError,
    #[error("io error")]
    IoError(#[from] std::io::Error),
    #[error("File encoding error")]
    EncodeError(#[from] std::string::FromUtf8Error),
}

pub trait FastaRecord {
    fn get_id(&self) -> &str;
    fn set_id(&mut self, id: String);
    fn id_mut(&mut self) -> &mut String;

    fn get_seq(&self) -> &str;
    fn set_seq(&mut self, seq: String);
    fn seq_mut(&mut self) -> &mut String;

    fn valid(&self) -> Result<(), FastaError>;
    fn empty(&self) -> bool;
    fn len(&self) -> usize;
    fn clear(&mut self);
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Record {
    id: String,
    seq: String,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct CowRecord<'a> {
    id: Cow<'a, String>,
    seq: Cow<'a, String>,
}

impl Record {
    pub fn new() -> Self {
        Record {
            id: String::from(""),
            seq: String::from(""),
        }
    }
}

impl FastaRecord for Record {
    fn valid(&self) -> Result<(), FastaError> {
        if !self.id.is_ascii() && self.seq.is_ascii() {
            return Err(FastaError::ValidationError("Non-ASCII character in record"));
        }

        if self.id.is_empty() {
            return Err(FastaError::ValidationError("id field cannot be empty"));
        }

        Ok(())
    }

    fn get_id(&self) -> &str {
        self.id.as_ref()
    }

    fn set_id(&mut self, id: String) {
        self.id = id;
    }

    fn id_mut(&mut self) -> &mut String {
        &mut self.id
    }

    fn get_seq(&self) -> &str {
        self.seq.as_ref()
    }

    fn set_seq(&mut self, seq: String) {
        self.seq = seq;
    }

    fn seq_mut(&mut self) -> &mut String {
        &mut self.seq
    }

    fn empty(&self) -> bool {
        self.id.is_empty() && self.seq.is_empty()
    }

    fn len(&self) -> usize {
        self.seq.len()
    }

    fn clear(&mut self) {
        self.id.clear();
        self.seq.clear();
    }
}

impl<'a> FastaRecord for CowRecord<'a> {
    fn valid(&self) -> Result<(), FastaError> {
        if !self.id.is_ascii() && self.seq.is_ascii() {
            return Err(FastaError::ValidationError("Non-ASCII character in record"));
        }

        if self.id.is_empty() {
            return Err(FastaError::ValidationError("id field cannot be empty"));
        }

        Ok(())
    }

    fn get_id(&self) -> &str {
        self.id.as_ref()
    }

    fn set_id(&mut self, id: String) {
        self.id = Cow::Owned(id);
    }

    fn id_mut(&mut self) -> &mut String {
        self.id.to_mut()
    }

    fn get_seq(&self) -> &str {
        self.seq.as_ref()
    }

    fn set_seq(&mut self, seq: String) {
        self.seq = Cow::Owned(seq);
    }

    fn seq_mut(&mut self) -> &mut String {
        self.seq.to_mut()
    }

    fn empty(&self) -> bool {
        self.id.is_empty() && self.seq.is_empty()
    }

    fn len(&self) -> usize {
        self.seq.len()
    }

    fn clear(&mut self) {
        self.id.to_mut().clear();
        self.seq.to_mut().clear();
    }
}

pub trait FastaReader {
    fn buffer(&mut self) -> &mut String;
    fn read_line(&mut self) -> Result<usize, std::io::Error>;
    fn inner(&mut self) -> &mut dyn BufRead;
}

#[derive(Debug)]
pub struct CheckedReader<T, R>
where
    T: BufRead,
    R: FastaRecord,
{
    inner: T,
    buffer: String,
    record_type: PhantomData<R>,
}

impl<T, R> FastaReader for CheckedReader<T, R>
where
    T: BufRead,
    R: FastaRecord,
{
    fn inner(&mut self) -> &mut dyn BufRead {
        &mut self.inner
    }

    fn buffer(&mut self) -> &mut String {
        &mut self.buffer
    }

    fn read_line(&mut self) -> Result<usize, std::io::Error> {
        self.inner.read_line(&mut self.buffer)
    }
}

impl<T, R> CheckedReader<T, R>
where
    T: BufRead,
    R: FastaRecord,
{
    pub fn new(f: T) -> Self {
        CheckedReader {
            inner: f,
            buffer: String::from(""),
            record_type: PhantomData,
        }
    }

    pub fn read_record(&mut self, record: &mut R) -> Result<(), FastaError> {
        record.clear();
        self.buffer().clear();

        if self.buffer.is_empty() {
            match self.inner.read_line(&mut self.buffer) {
                Ok(0) if record.empty() => {
                    return Ok(()); // EOF
                }
                Ok(_) => self.buffer.trim_end().to_string(),
                Err(e) => return Err(FastaError::IoError(e)),
            };
        }

        if !self.buffer.starts_with('>') {
            return Err(FastaError::MissingId);
        }

        record.set_id(self.buffer[1..].trim_end().to_string());
        self.buffer.clear();

        // read sequence content
        match self.inner.read_line(&mut self.buffer) {
            Ok(0) => return Err(FastaError::EofError),
            Ok(_) => (),
            Err(e) => return Err(FastaError::IoError(e)),
        };
        while !self.buffer.is_empty() && !self.buffer.starts_with('>') {
            record.seq_mut().push_str(&self.buffer);
            self.buffer.clear();
            match self.inner.read_line(&mut self.buffer) {
                Ok(_) => (),
                Err(e) => return Err(FastaError::IoError(e)),
            }
        }

        if record.get_seq().is_empty() {
            return Err(FastaError::MissingSequenceError);
        }

        record.seq_mut().retain(|c| !char::is_ascii_whitespace(&c));
        Ok(())
    }
}

impl<T, R> Iterator for CheckedReader<T, R>
where
    T: BufRead,
    R: FastaRecord + std::default::Default,
{
    type Item = Result<R, FastaError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Default::default();
        let res = self.read_record(&mut record);
        match res {
            Ok(()) if record.empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(e) => Some(Err(e)),
        }
    }
}

#[derive(Debug)]
pub struct Reader<T, R>
where
    T: BufRead,
    R: FastaRecord,
{
    inner: T,
    buffer: String,
    record_type: PhantomData<R>,
}

impl<T, R> Reader<T, R>
where
    T: BufRead,
    R: FastaRecord,
{
    pub fn new(f: T) -> Self {
        Reader {
            inner: f,
            buffer: String::from(""),
            record_type: PhantomData,
        }
    }

    pub fn inner(&mut self) -> &mut T {
        &mut self.inner
    }

    pub fn buffer(&mut self) -> &mut String {
        &mut self.buffer
    }

    pub fn read_line(&mut self) -> Result<usize, std::io::Error> {
        self.inner.read_line(&mut self.buffer)
    }

    pub fn read_record(&mut self, record: &mut R) {
        record.clear();

        if self.buffer.is_empty() {
            self.inner
                .read_line(&mut self.buffer)
                .expect("failed to read fasta line");
        }

        if !self.buffer.is_empty() {
            record.set_id(self.buffer[1..].trim_end().to_string());
        } else {
            // if there's no ID, we've reached EOF
            return;
        }
        self.buffer.clear();

        // read sequence content
        self.inner
            .read_line(&mut self.buffer)
            .expect("failed to read fasta line");
        while !self.buffer.is_empty() && !self.buffer.starts_with('>') {
            record.seq_mut().push_str(&self.buffer);
            self.buffer.clear();
            self.inner
                .read_line(&mut self.buffer)
                .expect("failed to read fasta line");
        }
        record.seq_mut().retain(|c| !char::is_ascii_whitespace(&c));
    }
}

impl<T, R> Iterator for Reader<T, R>
where
    T: BufRead,
    R: FastaRecord + std::default::Default,
{
    type Item = R;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Default::default();
        self.read_record(&mut record);
        if record.empty() {
            return None;
        }
        Some(record)
    }
}

#[cfg(test)]
mod tests {

    use crate::io::fasta::*;
    use crate::io::*;
    use std::fs::File;
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
        let mut reader = mmapped_fa_reader(&f).unwrap();
        let mut record = Record::new();
        for r in reader {
            assert!(r.valid().is_ok());
        }
    }

    #[test]
    #[should_panic]
    fn test_bad_fa_panics() {
        let fa_path = init_path("resources/test_data/corrupt.fa");
        let f = File::open(fa_path).unwrap();
        let reader = mmapped_fa_reader(&f).unwrap();
        for r in reader {
            eprintln!("{:?}", r);
        }
    }

    #[test]
    fn test_get_fields() {
        let fa_path = init_path("resources/test_data/test.fa");
        let f = File::open(fa_path).unwrap();
        let mut reader = mmapped_fa_reader(&f).unwrap();
        let record = reader.next().unwrap();

        assert!(record.get_id() == "SRR22092847.1.1");
        assert!(
            record.get_seq()
                == "GNTTAAAGCACATAAAGACAAATCGCTCCAGGGCAAAGNTTAAAGCACATAAAGACAAATCGCTCCAGGGCAAA"
        );
    }

    #[test]
    fn test_bad_fa_is_recoverable() {
        let fa_path = init_path("resources/test_data/trunc.fa");
        let f = File::open(fa_path).unwrap();
        let reader = mmapped_cfa_reader(&f).unwrap();
        for record in reader {
            match record {
                Ok(_) => {}
                Err(e) => eprintln!("{}", e),
            }
        }
    }
}
