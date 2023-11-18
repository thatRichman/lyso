use std::fmt::Display;
use std::str::Utf8Error;
use thiserror::Error;

pub(crate) mod parser;
pub mod reader;
// pub mod indexer;

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
    #[error("Error parsing fastq record")]
    ParseError,
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

impl Display for Record {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "@{} {}\n", self.id, self.desc)?;
        write!(f, "{}\n", self.seq)?;
        write!(f, "+\n")?;
        write!(f, "{}\n", self.qual)
    }
}
