#![feature(test)]

extern crate test;

#[cfg(test)]
mod benches {
    use lyso_fastq::{reader::FastqReader, FastqError, Record};
    use std::{fs::File, io::BufReader};

    use super::*;
    use test::{black_box, Bencher};

    #[bench]
    pub fn bench_read_fq(b: &mut Bencher) {
        let f = File::open("../benches/bench-fastq/med.fastq").unwrap();
        let mut reader = BufReader::new(&f);
        let mut fq_reader = black_box(FastqReader::new(reader));
        b.iter(|| {
            black_box(black_box(&mut fq_reader).collect::<Vec<Result<Record, FastqError>>>());
        });
    }
}
