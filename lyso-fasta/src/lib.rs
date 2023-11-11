use std::borrow::Cow;
use std::fmt::{write, Display};
use std::io::BufRead;
use std::iter::Iterator;
use thiserror::Error;

//pub mod indexer;
pub mod parser;
pub mod reader;

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
    #[error("Missing header")]
    TruncatedId,
    #[error("Parse error")]
    ParserError,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Record {
    id: String,
    seq: String,
}

impl<'a> Record {
    pub fn new() -> Self {
        Record {
            id: String::from(""),
            seq: String::from(""),
        }
    }
}

impl Display for Record {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\n", self.id)?;
        write!(f, "{}", self.seq)
    }
}
