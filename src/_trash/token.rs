fn readSeq(r: &mut impl BufRead) {

}

fn readQual(r: &mut impl BufRead) {

}

fn skipN(r: &mut impl BufRead, n: u64) {
    let mut _buf = String::from("");
    for i in 0..n {
        r.read_line(_buf);
    }
}


fn skip1(r: &mut impl BufRead) { 
        r.read_line(String::from(""))
}
