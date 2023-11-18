use nom::{
    bytes::complete::is_a as complete_is_a,
    bytes::streaming::{is_not, tag},
    combinator::{map_res, opt, cut},
    sequence::{pair, preceded, terminated, tuple},
    IResult,
};

#[inline]
fn start(input: &[u8]) -> IResult<&[u8], &[u8]> {
    tag(b"@")(input)
}

#[inline]
fn not_line_ending(input: &[u8]) -> IResult<&[u8], &[u8]> {
    is_not("\r\n")(input)
}

#[inline]
/// This uses the complete form of the `is_a` parser.
/// The reason for this is that streaming parsers make
/// it exceptionally difficult to differentiate between true EOF
/// and actually needing more data.
/// It is the responsibility of the reader implementing this parser
/// to ensure the passed buffer always ends on a newline.
fn line_ending(input: &[u8]) -> IResult<&[u8], &[u8]> {
    complete_is_a("\r\n")(input)
}

#[inline]
fn not_line_ending_or_space(input: &[u8]) -> IResult<&[u8], &[u8]> {
    is_not("\r\n ")(input)
}

#[inline]
fn header(input: &[u8]) -> IResult<&[u8], (&str, &str)> {
    let (i, (id, desc)) = terminated(
        pair(
            preceded(start, not_line_ending_or_space),
            opt(preceded(tag(" "), not_line_ending)),
        ),
        line_ending,
    )(input)?;
    Ok((
        i,
        (
            std::str::from_utf8(id).unwrap(),
            std::str::from_utf8(desc.unwrap_or(&[])).unwrap(),
        ),
    ))
}

#[inline]
fn line(input: &[u8]) -> IResult<&[u8], &str> {
    map_res(terminated(not_line_ending, line_ending), |x| {
        std::str::from_utf8(x)
    })(input)
}

#[inline]
fn comment(input: &[u8]) -> IResult<&[u8], &[u8]> {
    terminated(tag("+"), line)(input)
}

#[inline]
pub fn parse_record(input: &[u8]) -> IResult<&[u8], (&str, &str, &str, &str)> {
    let (i, ((id, desc), seq, _, qual)) = tuple((cut(header), line, comment, line))(input)?;
    Ok((i, (id, desc, seq, qual)))
}

#[cfg(test)]
mod tests {}
