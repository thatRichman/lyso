#![feature(test)]

extern crate test;

#[cfg(test)]
mod benches {
    use lyso_fasta::{reader::FastaReader, FastaError, Record};
    use std::{fs::File, io::BufReader};

    use super::*;
    use test::{black_box, Bencher};

    #[bench]
    pub fn bench_read_fa(b: &mut Bencher) {
        let f = File::open("../benches/bench-fasta/med.fa").unwrap();
        let mut reader = BufReader::new(&f);
        let mut fa_reader = FastaReader::new(reader).unwrap();
        b.iter(|| {
            black_box((&mut fa_reader).collect::<Vec<Result<Record, FastaError>>>());
        });
    }
}
