use fxhash::FxHashMap;
use nom::{
    bytes::complete::take_until,
    bytes::streaming::{tag, take},
    combinator::{map, map_parser},
    multi::{count, fill, length_data, many1},
    number::complete,
    number::streaming,
    sequence::{preceded, tuple},
    IResult,
};

use crate::{BamAuxField, BamAuxValue, BamHeader, BamReference, BamSeq, Record, BAM_MAGIC_STR};
use lyso_common::CigarOp;

// ============================== //
//    BEGIN BAM HEADER PARSING    //
// ============================== //

/// Parse BAM magic string
///
/// Attempts to match [66, 65, 77, 1].
pub fn bam_magic(input: &[u8]) -> IResult<&[u8], &[u8]> {
    tag(BAM_MAGIC_STR)(input)
}

/// Parse BAM header into tuple
///
/// Attempts to match BAM magic string (discarded), header text, and n_ref.
fn header(input: &[u8]) -> IResult<&[u8], (&[u8], u32)> {
    tuple((
        preceded(bam_magic, length_data(streaming::le_u32)),
        streaming::le_u32,
    ))(input)
}

/// Convert bytes into `BamHeader` struct
pub fn read_header(input: &[u8]) -> IResult<&[u8], BamHeader> {
    match header(input) {
        Ok((_i, (text_bytes, n_ref))) => IResult::Ok((
            _i,
            BamHeader {
                text: String::from_utf8_lossy(text_bytes).into_owned(),
                n_ref,
            },
        )),
        Err(e) => Err(e),
    }
}

// =========================== //
// BEGIN BAM REFERENCE PARSING //
// =========================== //

// nom streaming parsers are very unweildy for char pattern matching
// because they usually return Needed(1) so implement manually.
// TODO evaluate if this is faster than just using a regex

/// Validate reference name
///
/// SAMv1 1.2.1:
/// Reference name must be ASCII and match the following regex
/// [0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*
fn validate_ref_name(name: &str) -> Option<&str> {
    if name.chars().take(1).all(|c| !matches!(c, '=' | '*'))
        && name.chars().all(|c| {
            c.is_ascii_graphic()
                && !matches!(
                    c,
                    '\\' | '{' | '}' | '[' | ']' | '<' | '>' | '(' | ')' | ','
                )
        })
    {
        Some(name)
    } else {
        None
    }
}

fn reference_name(input: &[u8]) -> IResult<&[u8], &[u8]> {
    map_parser(length_data(streaming::le_u32), null_terminated_bytes)(input)
}

/// Parse BAM reference data into tuple
fn reference(input: &[u8]) -> IResult<&[u8], (&[u8], u32)> {
    tuple((reference_name, streaming::le_u32))(input)
}

/// Convert bytes into a BamReference
///
/// Attempts to parse `input` into BamReference, returning unconsumed input and BamReference if
/// successful.
pub fn read_reference(input: &[u8]) -> IResult<&[u8], BamReference> {
    match reference(input) {
        Ok((_i, (text_bytes, l_ref))) => IResult::Ok((
            _i,
            BamReference {
                name: validate_ref_name(std::str::from_utf8(text_bytes).unwrap())
                    .expect("Invalid reference name")
                    .to_string(),
                l_ref,
            },
        )),
        Err(e) => Err(e),
    }
}

/// Read `n` references into Vec<BamReference>
///
/// See also `read_reference`.
pub fn read_references<'a>(input: &'a [u8], buf: &mut [BamReference]) -> IResult<&'a [u8], ()> {
    fill(read_reference, buf)(input)
}

// =========================== //
//   BEGIN BAM CIGAR PARSING   //
// =========================== //

/// Converts unpacked CIGAR data into a single CigarOp
///
/// This expects a [u32; 2], as obtained by `unpack_cigar_op` parser.

fn to_cigar(input: [u32; 2]) -> CigarOp {
    match input[0] {
        0 => CigarOp::M(input[1]),
        1 => CigarOp::I(input[1]),
        2 => CigarOp::D(input[1]),
        3 => CigarOp::N(input[1]),
        4 => CigarOp::S(input[1]),
        5 => CigarOp::H(input[1]),
        6 => CigarOp::P(input[1]),
        7 => CigarOp::Eq(input[1]),
        8 => CigarOp::X(input[1]),
        otherwise => panic!("Invalid CigarOp {}{otherwise}", input[1]),
    }
}

/// Unpacks a compressed CIGAR operation
///
/// Reads a single u32 and unpacks operation + length.
/// See SAM v1 4.2

pub fn unpack_cigar_op(input: &[u8]) -> IResult<&[u8], [u32; 2]> {
    let (_i, v) = complete::le_u32(input)?;
    Ok((_i, [v & 4, v >> 4 | (v & 4)]))
}

/// Read bytes into vector of `CigarOp`s
///
/// Reads and unpacks `n_op` bytes, converting each to corresponding CigarOp variant.

pub fn read_cigar<'a>(input: &'a [u8], n_op: &u16) -> IResult<&'a [u8], Vec<CigarOp>> {
    let mut ops: Vec<CigarOp> = Vec::with_capacity(usize::try_from(*n_op).unwrap());
    let mut _i: &[u8] = input;
    for _ in 0..(*n_op) {
        (_i, _) = map(unpack_cigar_op, |v| ops.push(to_cigar(v)))(_i)?;
    }
    Ok((_i, ops))
}

// ============================== //
//   BEGIN BAM SEQUENCE PARSING   //
// ============================== //

/// Parse byte into BAM sequence
///
/// See SAM v1 4.2.3

pub fn to_sequence(input: &u8) -> BamSeq {
    match input {
        0 => BamSeq::Eq,
        1 => BamSeq::A,
        2 => BamSeq::C,
        3 => BamSeq::M,
        4 => BamSeq::G,
        5 => BamSeq::R,
        6 => BamSeq::S,
        7 => BamSeq::V,
        8 => BamSeq::T,
        9 => BamSeq::W,
        10 => BamSeq::Y,
        11 => BamSeq::H,
        12 => BamSeq::K,
        13 => BamSeq::D,
        14 => BamSeq::B,
        _ => BamSeq::N,
    }
}

/// Unpacks compressed BAM sequence
///
/// Each byte contains two sequence values.
/// Returns the new values as bytes.

fn unpack_sequence(input: &[u8]) -> IResult<&[u8], [u8; 2]> {
    let (_i, v) = complete::le_u8(input)?;
    Ok((_i, [v >> 4, v & 0x0F]))
}

/// Read bytes into vector of `Seq`s
///
/// The sequence field is bit-packed, thus the length of the returned vec is not `l_seq`
/// but rather (`l_seq` + 1) / 2. In the event that `l_seq` is odd, the final 4 bits are garbage
/// and automatically discarded.
pub fn read_sequence<'a>(input: &'a [u8], l_seq: &u32) -> IResult<&'a [u8], Vec<BamSeq>> {
    let mut seq: Vec<BamSeq> = Vec::with_capacity(usize::try_from((*l_seq + 1) / 2).unwrap());
    let mut _i: &[u8] = input;
    for _ in 0..seq.capacity() {
        (_i, _) = map(unpack_sequence, |v| {
            seq.push(to_sequence(&v[0]));
            seq.push(to_sequence(&v[1]));
        })(_i)?;
    }
    if l_seq % 2 != 0 {
        seq.pop();
    }
    Ok((_i, seq))
}

// ============================== //
//   BEGIN BAM SEQUENCE PARSING   //
// ============================== //

/// Read PHRED quality values
///
/// `n` is expected to be the value of BAM `seq_len` field.
fn read_quality(input: &[u8], n: u32) -> IResult<&[u8], Vec<u8>> {
    count(complete::le_u8, usize::try_from(n).unwrap())(input)
}

// ============================== //
//      BEGIN AUX BAM PARSING     //
// ============================== //

// TODO evaluate if there's really any benefit to returning [char; 2]
// instead of String.
/// Reads a two-character bam tag
fn bam_tag(input: &[u8]) -> IResult<&[u8], [char; 2]> {
    let mut buf: [u8; 2] = [0; 2];
    let (i, _) = fill(complete::le_u8, &mut buf)(input)?;
    Ok((i, [buf[0] as char, buf[1] as char]))
}

/// Parse bytes until encountering NULL (\0)
///
/// Consumes but does not return NULL.
fn null_terminated_bytes(input: &[u8]) -> IResult<&[u8], &[u8]> {
    let (mut i, r) = take_until(&[0u8] as &[u8])(input)?;
    (i, _) = take::<usize, &[u8], nom::error::Error<_>>(1usize)(i).unwrap();
    Ok((i, r))
}

/// Read a hex values into vector of u32s
///
/// Consumes all valid hex values.
fn hex_vec(input: &[u8]) -> IResult<&[u8], Vec<u32>> {
    many1(complete::hex_u32)(input)
}

/// Read variable-length auxilliary fields into BamAuxValue
///
/// Consumes subtype, length, and field, returning BamAuxValue.
fn aux_vec(input: &[u8]) -> IResult<&[u8], BamAuxValue> {
    let (i, (sub, len)) = tuple((complete::le_u8, complete::le_u32))(input)?;
    let len = usize::try_from(len).unwrap();
    match sub {
        b'c' => map(count(complete::le_i8, len), |v| BamAuxValue::Bc(v))(i),
        b'C' => map(count(complete::le_u8, len), |v| BamAuxValue::BC(v))(i),
        b's' => map(count(complete::le_i16, len), |v| BamAuxValue::Bs(v))(i),
        b'S' => map(count(complete::le_u16, len), |v| BamAuxValue::BS(v))(i),
        b'i' => map(count(complete::le_i32, len), |v| BamAuxValue::Bi(v))(i),
        b'I' => map(count(complete::le_u32, len), |v| BamAuxValue::BI(v))(i),
        b'f' => map(count(complete::le_f32, len), |v| BamAuxValue::Bf(v))(i),
        otherwise => panic!("Unknown BAM auxilliary field subtype {otherwise}"),
    }
}

/// Read BAM auxilliary fields into BamAuxField
///
/// Consumes tag, dtype, and value, returning BamAuxField
fn read_aux_field(input: &[u8]) -> IResult<&[u8], BamAuxField> {
    let (i, tag) = bam_tag(input)?;
    let (i, dtype) = complete::le_u8(i)?;
    let (i, value) = match dtype {
        b'A' => map(complete::le_u8, |v| BamAuxValue::from(v as char))(i)?,
        b'c' => map(complete::le_i8, |v| BamAuxValue::from(v))(i)?,
        b'C' => map(complete::le_u8, |v| BamAuxValue::from(v))(i)?,
        b's' => map(complete::le_i16, |v| BamAuxValue::from(v))(i)?,
        b'S' => map(complete::le_u16, |v| BamAuxValue::from(v))(i)?,
        b'i' => map(complete::le_i32, |v| BamAuxValue::from(v))(i)?,
        b'I' => map(complete::le_u32, |v| BamAuxValue::from(v))(i)?,
        b'f' => map(complete::le_f32, |v| BamAuxValue::from(v))(i)?,
        b'Z' => map(null_terminated_bytes, |v| {
            BamAuxValue::from(std::str::from_utf8(v).unwrap().to_owned())
        })(i)?,
        b'H' => map(hex_vec, |v| BamAuxValue::H(v))(i)?,
        b'B' => aux_vec(i)?,
        otherwise => panic!("Invalid BAM auxilliary field {otherwise}"),
    };
    Ok((i, BamAuxField { tag, value }))
}

// ============================== //
//   BEGIN BAM ALIGNMENT PARSING   //
// ============================== //

// The following parsers are strictly unnecessary
// but naming them makes me happy

/// parse block size
pub fn block_size(input: &[u8]) -> IResult<&[u8], u32> {
    complete::le_u32(input)
}

/// parse ref_id
fn ref_id(input: &[u8]) -> IResult<&[u8], i32> {
    complete::le_i32(input)
}

/// parse pos
fn pos(input: &[u8]) -> IResult<&[u8], i32> {
    complete::le_i32(input)
}

/// parse l_read_name
fn l_read_name(input: &[u8]) -> IResult<&[u8], u8> {
    complete::le_u8(input)
}

/// parse mapq
fn mapq(input: &[u8]) -> IResult<&[u8], u8> {
    complete::le_u8(input)
}

/// parse bin
fn bin(input: &[u8]) -> IResult<&[u8], u16> {
    complete::le_u16(input)
}

/// parse n_cigar_op
pub fn n_cigar_op(input: &[u8]) -> IResult<&[u8], u16> {
    complete::le_u16(input)
}

/// parse flag
fn flag(input: &[u8]) -> IResult<&[u8], u16> {
    complete::le_u16(input)
}

/// parse l_seq
fn l_seq(input: &[u8]) -> IResult<&[u8], u32> {
    complete::le_u32(input)
}

/// parse next_ref_id
fn next_ref_id(input: &[u8]) -> IResult<&[u8], i32> {
    complete::le_i32(input)
}

/// parse next_pos
fn next_pos(input: &[u8]) -> IResult<&[u8], i32> {
    complete::le_i32(input)
}

/// parse tlen
fn tlen(input: &[u8]) -> IResult<&[u8], i32> {
    complete::le_i32(input)
}

/// parse read_name
///
/// n is expected to be value parsed from `l_read_name`
fn read_name(input: &[u8], n: u8) -> IResult<&[u8], &[u8]> {
    take(n)(input)
}

/// Convert Vec<BamAuxField> to HashMap
///
/// Maps BamAuxField.tag (as String) to BamAuxField.

fn aux_to_hash(fields: Vec<BamAuxField>) -> FxHashMap<String, BamAuxField> {
    let mut hmap = FxHashMap::default();
    for f in fields {
        hmap.insert(f.tag.iter().collect(), f);
    }
    hmap
}

/// Maybe correct for long CIGAR fields
///
/// If the criteria described in SAMv1 4.2.2 are met,
/// update `n_cigar_op` and `cigar_op` fields, and remove the
/// "CG" aux field.

fn maybe_correct_cigar(
    n_cigar_op: &mut u16,
    seq_len: &usize,
    cigar: &mut Vec<CigarOp>,
    aux_hash: &mut FxHashMap<String, BamAuxField>,
    reference: &BamReference,
) {
    if *n_cigar_op == 2
        && aux_hash.contains_key("CG")
        && cigar
            == &[
                CigarOp::S(u32::try_from(*seq_len).unwrap()),
                CigarOp::N(reference.l_ref),
            ]
    {
        match aux_hash.get("CG") {
            Some(BamAuxField {
                tag: _,
                value: BamAuxValue::BI(v),
            }) => {
                *n_cigar_op = u16::try_from(v.len()).unwrap();
                *cigar = v
                    .chunks_exact(2)
                    .map(|v: &[u32]| to_cigar(v.try_into().unwrap()))
                    .collect::<Vec<CigarOp>>();
                aux_hash.remove("CG");
            }
            _ => {}
        }
    }
}

/// Read a complete alignment record
pub fn read_alignment<'a>(
    input: &'a [u8],
    references: &Vec<BamReference>,
) -> IResult<&'a [u8], Record> {
    let (
        i,
        (
            block_size,
            ref_id,
            pos,
            l_read_name,
            mapq,
            bin,
            mut n_cigar_op,
            flag,
            l_seq,
            next_ref_id,
            next_pos,
            tlen,
        ),
    ) = tuple((
        block_size,
        ref_id,
        pos,
        l_read_name,
        mapq,
        bin,
        n_cigar_op,
        flag,
        l_seq,
        next_ref_id,
        next_pos,
        tlen,
    ))(input)?;

    // each of these requires one of the above items
    let (i, read_name_bytes) = read_name(i, l_read_name)?;
    let read_name = String::from_utf8_lossy(read_name_bytes).to_string();
    let (i, mut cigar) = read_cigar(i, &n_cigar_op)?;
    let (i, seq) = read_sequence(i, &l_seq)?;

    let (mut i, raw_qual) = read_quality(i, l_seq)?;
    let qual = if raw_qual.iter().all(|q| q == &255) {
        None
    } else {
        Some(raw_qual)
    };

    let mut aux_fields: Vec<BamAuxField> = Vec::with_capacity(4);
    if !i.is_empty() {
        (i, aux_fields) = many1(read_aux_field)(i).unwrap();
    }
    let mut aux_hash: Option<FxHashMap<String, BamAuxField>> = if aux_fields.len() > 0 {
        Some(aux_to_hash(aux_fields))
    } else {
        None
    };

    maybe_correct_cigar(
        &mut n_cigar_op,
        &seq.len(),
        &mut cigar,
        aux_hash.as_mut().unwrap(),
        &references[usize::try_from(ref_id).unwrap()],
    );

    Ok((
        i,
        Record {
            block_size,
            ref_id,
            pos,
            l_read_name,
            mapq,
            bin,
            n_cigar_op,
            flag,
            l_seq,
            next_ref_id,
            next_pos,
            tlen,
            read_name,
            cigar,
            seq,
            qual,
            aux: aux_hash,
        },
    ))
}
