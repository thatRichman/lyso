pub mod indexer;
pub mod parser;
pub mod reader;

use fxhash::FxHashMap;
use lyso_common::CigarOp;
use nom::error::ParseError;
use std::fmt::{self, Display};
use thiserror::Error;

const BAM_MAGIC_STR: [u8; 4] = [66, 65, 77, 1];
const MAX_BLOCK_SIZE: usize = 65536;

#[derive(Debug)]
/// Sequence primitives
/// See SAM v1 section 4.2
pub enum BamSeq {
    Eq,
    A,
    C,
    M,
    G,
    R,
    S,
    V,
    T,
    W,
    Y,
    H,
    K,
    D,
    B,
    N,
}

impl Display for BamSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BamSeq::Eq => write!(f, "="),
            BamSeq::A => write!(f, "A"),
            BamSeq::C => write!(f, "C"),
            BamSeq::M => write!(f, "M"),
            BamSeq::N => write!(f, "N"),
            BamSeq::G => write!(f, "G"),
            BamSeq::R => write!(f, "R"),
            BamSeq::V => write!(f, "V"),
            BamSeq::T => write!(f, "T"),
            BamSeq::B => write!(f, "B"),
            BamSeq::W => write!(f, "W"),
            BamSeq::Y => write!(f, "Y"),
            BamSeq::S => write!(f, "S"),
            BamSeq::K => write!(f, "K"),
            BamSeq::H => write!(f, "H"),
            BamSeq::D => write!(f, "D"),
        }
    }
}

#[derive(Error, Debug)]
pub enum BamError {
    #[error("Unexpected EOF")]
    EofError,
    #[error("Missing BAM Magic String")]
    MissingMagicString,
    #[error("I/O error")]
    IoError(#[from] std::io::Error),
    #[error("File encoding error")]
    EncodeError(#[from] std::string::FromUtf8Error),
    #[error("Parse error")]
    ParseError,
    #[error("TryFromInt Error")]
    TryFromInt(#[from] std::num::TryFromIntError),
}

/// Auxilliary BAM field
///
/// arbitrary tag names are supported but must be of length 2
/// See BamAuxValue for possible value types.
#[derive(Debug)]
pub struct BamAuxField {
    tag: [char; 2],
    value: BamAuxValue,
}

impl Display for BamAuxField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}:{}", self.tag[0], self.tag[1], self.value)
    }
}

/// Auxilliary BAM value encodings
///
/// Display implementation will write in SAM format.
/// See SAM v1 section 4.2.4
#[allow(non_camel_case_types)]
#[derive(Debug)]
pub enum BamAuxValue {
    A(char),
    c(i8),
    C(u8),
    s(i16),
    S(u16),
    i(i32),
    I(u32),
    f(f32),
    Z(String),
    H(Vec<u32>),
    Bc(Vec<i8>),
    BC(Vec<u8>),
    Bs(Vec<i16>),
    BS(Vec<u16>),
    Bi(Vec<i32>),
    BI(Vec<u32>),
    Bf(Vec<f32>),
}

/// All integer types are 'i' in SAM format
impl Display for BamAuxValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BamAuxValue::A(v) => write!(f, "A:{v}"),
            BamAuxValue::c(v) => write!(f, "i:{v}"),
            BamAuxValue::C(v) => write!(f, "i:{v}"),
            BamAuxValue::s(v) => write!(f, "i:{v}"),
            BamAuxValue::S(v) => write!(f, "i:{v}"),
            BamAuxValue::i(v) => write!(f, "i:{v}"),
            BamAuxValue::I(v) => write!(f, "i:{v}"),
            BamAuxValue::Z(v) => write!(f, "Z:{v}"),
            _ => todo!(),
        }
    }
}

impl From<u32> for BamAuxValue {
    fn from(value: u32) -> Self {
        BamAuxValue::I(value)
    }
}

impl From<i32> for BamAuxValue {
    fn from(value: i32) -> Self {
        BamAuxValue::i(value)
    }
}

impl From<u8> for BamAuxValue {
    fn from(value: u8) -> Self {
        BamAuxValue::C(value)
    }
}

impl From<i8> for BamAuxValue {
    fn from(value: i8) -> Self {
        BamAuxValue::c(value)
    }
}

impl From<u16> for BamAuxValue {
    fn from(value: u16) -> Self {
        BamAuxValue::S(value)
    }
}

impl From<i16> for BamAuxValue {
    fn from(value: i16) -> Self {
        BamAuxValue::s(value)
    }
}

impl From<f32> for BamAuxValue {
    fn from(value: f32) -> Self {
        BamAuxValue::f(value)
    }
}

impl From<char> for BamAuxValue {
    fn from(value: char) -> Self {
        BamAuxValue::A(value)
    }
}

impl From<String> for BamAuxValue {
    fn from(value: String) -> Self {
        BamAuxValue::Z(value)
    }
}

impl From<Vec<u8>> for BamAuxValue {
    fn from(value: Vec<u8>) -> Self {
        BamAuxValue::BC(value)
    }
}

impl From<Vec<i8>> for BamAuxValue {
    fn from(value: Vec<i8>) -> Self {
        BamAuxValue::Bc(value)
    }
}

impl From<Vec<u16>> for BamAuxValue {
    fn from(value: Vec<u16>) -> Self {
        BamAuxValue::BS(value)
    }
}

impl From<Vec<i16>> for BamAuxValue {
    fn from(value: Vec<i16>) -> Self {
        BamAuxValue::Bs(value)
    }
}

impl From<Vec<u32>> for BamAuxValue {
    fn from(value: Vec<u32>) -> Self {
        BamAuxValue::BI(value)
    }
}

impl From<Vec<i32>> for BamAuxValue {
    fn from(value: Vec<i32>) -> Self {
        BamAuxValue::Bi(value)
    }
}

impl From<Vec<f32>> for BamAuxValue {
    fn from(value: Vec<f32>) -> Self {
        BamAuxValue::Bf(value)
    }
}

/// A BAM alignment record
#[derive(Debug, Default)]
pub struct Record {
    block_size: u32,
    ref_id: i32,
    ref_name: String,
    pos: i32,
    l_read_name: u8,
    mapq: u8,
    bin: u16,
    n_cigar_op: u16,
    flag: u16,
    l_seq: u32,
    next_ref_id: i32,
    next_ref_name: String,
    next_pos: i32,
    tlen: i32,
    read_name: String,
    cigar: Vec<CigarOp>,
    seq: Vec<BamSeq>,
    qual: Option<Vec<u8>>,
    aux: Option<FxHashMap<String, BamAuxField>>, // everything else
}

impl Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.read_name,
            self.flag,
            self.ref_name,
            self.pos + 1, // SAM is 1-based
            self.mapq,
            self.cigar.iter().map(|x| x.to_string()).collect::<String>(),
            self.next_ref_name,
            self.next_pos + 1,
            self.tlen,
            self.seq.iter().map(|x| x.to_string()).collect::<String>(),
            std::str::from_utf8(
                self.qual
                    .as_ref()
                    .unwrap_or(&vec![42u8; 1])
                    .iter()
                    .map(|x| x + 33)
                    .collect::<Vec<u8>>()
                    .as_ref()
            )
            .unwrap_or("*")
        )
        .unwrap();
        if self.aux.is_some() {
            for val in self.aux.as_ref().unwrap().values() {
                write!(f, "\t{val}").unwrap();
            }
        }
        Ok(())
    }
}

/// Representation of BAM Reference record
///
/// Display implementation will write in SAM format.
#[derive(Debug)]
pub struct BamReference {
    name: String,
    l_ref: u32,
}

/// Representation of BAM header field
///
/// Display implementation will write in SAM format.
#[derive(Debug)]
pub struct BamHeader {
    text: String,
    n_ref: u32,
}

#[derive(Clone, Copy, Debug, Default)]
pub enum PhredEncoding {
    #[default]
    Phred33 = 33,
    Phred64 = 64,
    Unknown = 0,
}

pub fn guess_phred_encoding(scores: &[u8]) -> PhredEncoding {
    let min = scores.iter().min().unwrap_or(&0);
    let max = scores.iter().max().unwrap_or(&0);
    if min < &59 && max <= &74 {
        return PhredEncoding::Phred33;
    }
    if min >= &64 && max > &73 {
        return PhredEncoding::Phred64;
    }
    PhredEncoding::Unknown
}
