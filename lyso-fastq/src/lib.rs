use std::borrow::Cow;
use std::io::{BufRead, Seek, SeekFrom};
use std::iter::Iterator;
use std::marker::PhantomData;
use std::str::Utf8Error;
use thiserror::Error;

pub mod index;

#[derive(Error, Debug)]
pub enum FastqError {
    #[error("fastq validation error")]
    ValidationError(&'static str),
    #[error("end of file error")]
    EofError,
    #[error("missing id error")]
    MissingId,
    #[error("truncated id error")]
    TruncatedId,
    #[error("sequence-quality length mismatch")]
    SeqQualMismatch,
    #[error("index mismatch error")]
    IndexMismatch,
    #[error("io error")]
    IoError(#[from] std::io::Error),
    #[error("file encoding error")]
    EncodeError(#[from] Utf8Error),
}

/// The FastqRecord provides a standard contract
/// for interfacing with any struct that implements FastqRecord
pub trait FastqRecord {
    fn get_id(&self) -> &str;
    fn id_mut(&mut self) -> &mut String;
    fn get_desc(&self) -> &str;

    fn get_seq(&self) -> &str;
    fn seq_mut(&mut self) -> &mut String;

    fn get_qual(&self) -> &str;
    fn qual_mut(&mut self) -> &mut String;

    fn set_id(&mut self, id: String);
    fn set_seq(&mut self, id: String);
    fn set_qual(&mut self, id: String);
    fn set_desc(&mut self, id: String);

    fn valid(&self) -> Result<(), FastqError>;
    fn empty(&self) -> bool;
    fn len(&self) -> usize;
    fn clear(&mut self);
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Record {
    id: String,
    desc: String,
    seq: String,
    qual: String,
}

impl Record {
    pub fn new() -> Self {
        Record {
            id: String::from(""),
            desc: String::from(""),
            seq: String::from(""),
            qual: String::from(""),
        }
    }
}

impl FastqRecord for Record {
    fn valid(&self) -> Result<(), FastqError> {
        if !self.id.is_ascii()
            && self.desc.is_ascii()
            && self.seq.is_ascii()
            && self.qual.is_ascii()
            && self.seq.is_ascii()
        {
            return Err(FastqError::ValidationError("Non-ASCII character in record"));
        }

        if self.id.is_empty() {
            return Err(FastqError::ValidationError("id field cannot be empty"));
        }

        if !self.seq.len() == self.qual.len() {
            return Err(FastqError::SeqQualMismatch);
        }

        Ok(())
    }

    fn get_id(&self) -> &str {
        self.id.as_ref()
    }

    fn id_mut(&mut self) -> &mut String {
        &mut self.id
    }

    fn get_seq(&self) -> &str {
        self.seq.as_ref()
    }

    fn seq_mut(&mut self) -> &mut String {
        &mut self.seq
    }

    fn get_qual(&self) -> &str {
        self.qual.as_ref()
    }

    fn qual_mut(&mut self) -> &mut String {
        &mut self.qual
    }

    fn get_desc(&self) -> &str {
        self.desc.as_ref()
    }

    fn set_id(&mut self, id: String) {
        self.id = id;
    }

    fn set_seq(&mut self, seq: String) {
        self.seq = seq;
    }

    fn set_qual(&mut self, qual: String) {
        self.qual = qual;
    }

    fn set_desc(&mut self, desc: String) {
        self.desc = desc;
    }

    fn empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_empty() && self.seq.is_empty() && self.qual.is_empty()
    }

    fn len(&self) -> usize {
        self.seq.len()
    }

    fn clear(&mut self) {
        self.id.clear();
        self.desc.clear();
        self.seq.clear();
        self.qual.clear();
    }
}

/// CowRecords store fastq content as Cow<String>
/// rather than String. This type of record is
/// best used when you are processing a large
/// number of reads but will only rarely be editing
/// fields in those reads. You can alternatively
/// use Cow<Record> in your functions, but this means
/// the entire record will be cloned when you modify
/// any field in the record.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CowRecord<'a> {
    id: Cow<'a, String>,
    desc: Cow<'a, String>,
    seq: Cow<'a, String>,
    qual: Cow<'a, String>,
}

impl<'a> CowRecord<'a> {
    pub fn new() -> Self {
        CowRecord {
            id: Cow::Owned(String::from("")),
            desc: Cow::Owned(String::from("")),
            seq: Cow::Owned(String::from("")),
            qual: Cow::Owned(String::from("")),
        }
    }
}

impl<'a> FastqRecord for CowRecord<'a> {
    fn valid(&self) -> Result<(), FastqError> {
        if !self.id.is_ascii()
            && self.desc.is_ascii()
            && self.seq.is_ascii()
            && self.qual.is_ascii()
            && self.seq.is_ascii()
        {
            return Err(FastqError::ValidationError("Non-ASCII character in record"));
        }

        if self.id.is_empty() {
            return Err(FastqError::ValidationError("id field cannot be empty"));
        }

        if !self.seq.len() == self.qual.len() {
            return Err(FastqError::SeqQualMismatch);
        }

        Ok(())
    }

    fn get_id(&self) -> &str {
        self.id.as_ref()
    }

    fn id_mut(&mut self) -> &mut String {
        self.id.to_mut()
    }

    fn get_seq(&self) -> &str {
        self.seq.as_ref()
    }

    fn seq_mut(&mut self) -> &mut String {
        self.seq.to_mut()
    }

    fn get_qual(&self) -> &str {
        self.qual.as_ref()
    }

    fn qual_mut(&mut self) -> &mut String {
        self.qual.to_mut()
    }

    fn get_desc(&self) -> &str {
        self.desc.as_ref()
    }

    fn set_id(&mut self, id: String) {
        self.id = Cow::Owned(id);
    }

    fn set_seq(&mut self, seq: String) {
        self.seq = Cow::Owned(seq);
    }

    fn set_qual(&mut self, qual: String) {
        self.qual = Cow::Owned(qual);
    }

    fn set_desc(&mut self, desc: String) {
        self.desc = Cow::Owned(desc);
    }

    fn empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_empty() && self.seq.is_empty() && self.qual.is_empty()
    }

    fn len(&self) -> usize {
        self.seq.len()
    }

    fn clear(&mut self) {
        self.id.to_mut().clear();
        self.desc.to_mut().clear();
        self.seq.to_mut().clear();
        self.qual.to_mut().clear();
    }
}

pub trait FastqReader {
    fn buffer(&mut self) -> &mut String;
    fn read_line(&mut self) -> Result<usize, std::io::Error>;
    fn inner(&mut self) -> &mut dyn BufRead;
    fn set_position(&mut self, pos: u64) -> Result<u64, std::io::Error>;
    fn position(&mut self) -> Result<u64, std::io::Error>;
}

/// This is a checked fastq reader.
/// This variant will return Result<Option<Record>, FastqError>,
/// where FastqError is an enum of possible error conditions that
/// can occur when parsing a fastq file. This form of reader should
/// only be used when it is absolutely necessary to be able to recover
/// from a malformed input, because there is a performance hit
/// incurred from validation.
pub struct CheckedReader<T, R>
where
    T: BufRead,
    R: FastqRecord,
{
    inner: T,
    buffer: String,
    record_type: PhantomData<R>,
}

impl<T, R> FastqReader for CheckedReader<T, R>
where
    T: BufRead + Seek,
    R: FastqRecord,
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

    fn set_position(&mut self, pos: u64) -> Result<u64, std::io::Error> {
        self.inner.seek(SeekFrom::Start(pos))
    }

    fn position(&mut self) -> Result<u64, std::io::Error> {
        self.inner.stream_position()
    }
}

impl<T, R> CheckedReader<T, R>
where
    T: BufRead,
    R: FastqRecord,
{
    pub fn new(f: T) -> Self {
        CheckedReader {
            inner: f,
            buffer: String::from(""),
            record_type: PhantomData,
        }
    }

    pub fn read_record(&mut self, record: &mut R) -> Result<(), FastqError> {
        record.clear();
        self.buffer.clear();

        match self.inner.read_line(&mut self.buffer) {
            Ok(0) if record.empty() => return Ok(()), // EOF
            Ok(_) => {}
            Err(e) => return Err(FastqError::IoError(e)),
        };

        if !self.buffer.starts_with('@') {
            return Err(FastqError::MissingId);
        }

        // assume all content after first whitespace is description
        let mut header = self.buffer[1..].trim_end().splitn(2, char::is_whitespace);
        match header.next() {
            Some(v) => record.set_id(v.to_string()),
            None => return Err(FastqError::TruncatedId),
        }
        record.set_desc(header.next().unwrap_or("").to_string());
        self.buffer.clear();

        // read sequence content, track # of lines it takes up
        self.inner.read_line(&mut self.buffer)?;
        let mut seq_lines = 0;
        while !self.buffer.is_empty() && !self.buffer.starts_with('+') {
            record.seq_mut().push_str(self.buffer.trim_end());
            self.buffer.clear();
            match self.inner.read_line(&mut self.buffer) {
                Ok(0) => return Err(FastqError::EofError),
                Ok(_) => (),
                Err(e) => return Err(FastqError::IoError(e)),
            }
            seq_lines += 1;
        }

        // read the same # of lines for quality
        self.buffer.clear();
        for _ in 0..seq_lines {
            match self.inner.read_line(&mut self.buffer) {
                Ok(0) => return Err(FastqError::EofError),
                Err(e) => return Err(FastqError::IoError(e)),
                Ok(_) => {}
            }
        }
        self.buffer.retain(|c| !char::is_ascii_whitespace(&c));
        record.set_qual(self.buffer.clone());
        if record.get_qual().len() != record.get_seq().len() {
            return Err(FastqError::SeqQualMismatch);
        }

        Ok(())
    }
}

impl<T, R> Iterator for CheckedReader<T, R>
where
    T: BufRead,
    R: FastqRecord + std::default::Default,
{
    type Item = Result<R, FastqError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Default::default();
        match self.read_record(&mut record) {
            Ok(()) if record.empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(e) => Some(Err(e)),
        }
    }
}

pub struct Reader<T, R>
where
    T: BufRead,
    R: FastqRecord,
{
    inner: T,
    buffer: String,
    record_type: PhantomData<R>,
}

/// This is an unchecked fastq reader. This means it does not perform
/// error handling. If the reader encounters an error while parsing
/// a fastq record, it will panic. No gaurantee is made that returned
/// records will be complete or valid. If desired, you can still
/// call valid() on individual records produced by this reader.

impl<T, R> FastqReader for Reader<T, R>
where
    T: BufRead + Seek,
    R: FastqRecord,
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

    fn set_position(&mut self, pos: u64) -> Result<u64, std::io::Error> {
        self.inner.seek(SeekFrom::Start(pos))
    }

    fn position(&mut self) -> Result<u64, std::io::Error> {
        self.inner.stream_position()
    }
}

impl<T, R> Reader<T, R>
where
    T: BufRead,
    R: FastqRecord,
{
    pub fn new(f: T) -> Self {
        Reader {
            inner: f,
            buffer: String::from(""),
            record_type: PhantomData,
        }
    }

    pub fn read_record(&mut self, record: &mut R) {
        record.clear();
        self.buffer.clear();

        self.inner
            .read_line(&mut self.buffer)
            .expect("failed to read fastq line");

        if self.buffer.is_empty() {
            return;
        }

        // assume all content after first whitespace is description
        let mut header = self.buffer[1..].trim_end().splitn(2, char::is_whitespace);
        record.set_id(header.next().expect("").to_string());
        record.set_desc(header.next().expect("").to_string());
        self.buffer.clear();

        // read sequence content, track # of lines it takes up
        self.inner.read_line(&mut self.buffer).expect("");
        let mut seq_lines = 0;
        while !self.buffer.is_empty() && !self.buffer.starts_with('+') {
            record.seq_mut().push_str(self.buffer.trim_end());
            self.buffer.clear();
            self.inner
                .read_line(&mut self.buffer)
                .expect("failed to read fastq line");
            seq_lines += 1;
        }

        // read the same # of lines for quality
        self.buffer.clear();
        for _ in 0..seq_lines {
            self.inner
                .read_line(&mut self.buffer)
                .expect("failed to read fastq line");
        }
        self.buffer.retain(|c| !char::is_ascii_whitespace(&c));
        record.set_qual(self.buffer.clone());
    }
}

impl<T, R> Iterator for Reader<T, R>
where
    T: BufRead,
    R: FastqRecord + std::default::Default,
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

    use crate::*;
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
        let reader: Reader<_, Record> = Reader::new(b);
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
        let reader: Reader<_, Record> = Reader::new(b);
        for _ in reader {
            continue;
        }
    }

    #[test]
    fn test_get_fields() {
        let fq_path = init_path("resources/test_data/test.fastq");
        let f = File::open(fq_path).unwrap();
        let b = BufReader::new(f);
        let mut reader: Reader<_, Record> = Reader::new(b);
        let record = reader.next().expect("bad record!");

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
        let mut reader: CheckedReader<_, Record> = CheckedReader::new(b);
        assert!(reader.all(|x| x.expect("failed to read record").valid().is_ok()))
    }

    #[test]
    fn test_bad_fq_is_recoverable() {
        let fq_path = init_path("resources/test_data/trunc.fastq");
        let f = File::open(fq_path).unwrap();
        let b = BufReader::new(f);
        let reader: CheckedReader<_, Record> = CheckedReader::new(b);
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
        let reader: CheckedReader<_, Record> = CheckedReader::new(b);
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
