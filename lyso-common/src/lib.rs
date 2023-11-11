use std::fmt::{self, Display};

pub mod util;

#[derive(Debug, PartialEq)]
// CIGAR operations
// See SAM v1 section 1.4.6
pub enum CigarOp {
    M(u32),
    I(u32),
    D(u32),
    N(u32),
    S(u32),
    H(u32),
    P(u32),
    Eq(u32),
    X(u32),
}

impl Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CigarOp::M(v) => write!(f, "{v}M"),
            CigarOp::I(v) => write!(f, "{v}I"),
            CigarOp::D(v) => write!(f, "{v}D"),
            CigarOp::N(v) => write!(f, "{v}N"),
            CigarOp::S(v) => write!(f, "{v}S"),
            CigarOp::H(v) => write!(f, "{v}H"),
            CigarOp::P(v) => write!(f, "{v}P"),
            CigarOp::Eq(v) => write!(f, "{v}="),
            CigarOp::X(v) => write!(f, "{v}X"),
        }
    }
}
