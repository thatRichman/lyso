use nom::{
    bytes::streaming::{is_a, is_not, tag},
    combinator::{map, map_res},
    sequence::{pair, preceded, terminated},
    IResult,
};

fn start(input: &[u8]) -> IResult<&[u8], &[u8]> {
    tag(b">")(input)
}

fn not_line_ending(input: &[u8]) -> IResult<&[u8], &[u8]> {
    is_not("\r\n")(input)
}

fn line_ending(input: &[u8]) -> IResult<&[u8], &[u8]> {
    is_a("\r\n")(input)
}

fn header(input: &[u8]) -> IResult<&[u8], String> {
    map_res(
        terminated(preceded(start, not_line_ending), line_ending),
        |x| String::from_utf8(x.to_vec()),
    )(input)
}

fn seq(input: &[u8]) -> IResult<&[u8], &[u8]> {
    is_not(">")(input)
}

fn merge_seq(s: &str) -> String {
    s.replace("\n", "")
}

fn merge_seq_2(s: &str) -> String {
    s.split("\n").collect::<String>()
}

fn sequence(input: &[u8]) -> IResult<&[u8], String> {
    map(seq, |x| merge_seq_2(std::str::from_utf8(x).unwrap()))(input)
}

pub fn parse_record(input: &[u8]) -> IResult<&[u8], (String, String)> {
    pair(header, sequence)(input)
}
