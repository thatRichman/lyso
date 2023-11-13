use nom::{
    bytes::complete::is_not,
    bytes::streaming::{is_a, tag},
    combinator::{map, map_res},
    sequence::{pair, preceded, terminated},
    IResult,
};

#[inline]
fn start(input: &[u8]) -> IResult<&[u8], &[u8]> {
    tag(b">")(input)
}

#[inline]
fn not_line_ending(input: &[u8]) -> IResult<&[u8], &[u8]> {
    is_not("\r\n")(input)
}

#[inline]
fn line_ending(input: &[u8]) -> IResult<&[u8], &[u8]> {
    is_a("\r\n")(input)
}

#[inline]
fn header(input: &[u8]) -> IResult<&[u8], String> {
    map_res(
        terminated(preceded(start, not_line_ending), line_ending),
        |x| String::from_utf8(x.to_vec()),
    )(input)
}

/// !IMPORTANT!
/// This parser uses bytes::complete::is_not.
/// Thus, you must be certain that you have read the entire
/// sequence content into the input before parsing or else
/// the sequence will be truncated. The reason for doing this
/// is because nom streaming parsers make it exceptionally
/// inconvenient to handle the distinction between actual EOF
/// and just needing more data. It is the responisibility
/// of the reader implementing these parsers to handle this.
/// `FastaReader` in this library does so by always reading to the
/// next '>' before parsing.
#[inline]
fn seq(input: &[u8]) -> IResult<&[u8], &[u8]> {
    is_not(">")(input)
}

#[inline]
fn remove_newlines(s: &str) -> String {
    s.split(['\r', '\n']).collect::<String>()
}

#[inline]
fn sequence(input: &[u8]) -> IResult<&[u8], String> {
    map(seq, |x| remove_newlines(std::str::from_utf8(x).unwrap()))(input)
}

#[inline]
pub fn parse_record(input: &[u8]) -> IResult<&[u8], (String, String)> {
    pair(header, sequence)(input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_start() {
        assert!(start(&[b'>',]).is_ok())
    }

    #[test]
    fn test_not_line_ending() {
        assert!(
            not_line_ending(b"ABCZzi8535!)&^%$!@$~&^\n") == Ok((b"\n", b"ABCZzi8535!)&^%$!@$~&^"))
        );
        assert!(
            not_line_ending(b"ABCZzi8535!)&^%$!@$~&^\r") == Ok((b"\r", b"ABCZzi8535!)&^%$!@$~&^"))
        );
    }

    #[test]
    fn test_line_ending() {
        assert!(line_ending(b"\n ") == Ok((b" ", b"\n")));
        assert!(line_ending(b"\r\n ") == Ok((b" ", b"\r\n")));
    }

    #[test]
    fn test_header() {
        assert!(header(b">SRR 123\n ") == Ok((b" ", String::from("SRR 123"))))
    }

    #[test]
    fn test_seq() {
        assert!(seq(b"ATGCN\n>") == Ok((b">", b"ATGCN\n")))
    }

    #[test]
    fn test_sequence() {
        assert!(sequence(b"ATGCN\nATGCN") == Ok((&[], String::from("ATGCNATGCN"))))
    }

    #[test]
    fn test_parse_record() {
        assert!(
            parse_record(b">A\nATGCN\n") == Ok((&[], (String::from("A"), String::from("ATGCN"))))
        );
        assert!(
            parse_record(b">B\nATGCN") == Ok((&[], (String::from("B"), String::from("ATGCN"))))
        );
    }
}
