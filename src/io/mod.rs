pub mod fasta;
pub mod fastq;

use std::fs::File;
use std::io::Cursor;
use std::path::Path;

use memmap2::{Advice, Mmap};

const GZB1: u8 = 0x1f;
const GZB2: u8 = 0x8b;

/// UNSAFE
/// given a file path, create a memory map to file contents
pub fn open_mmapped<P: AsRef<Path>>(fpath: P) -> Result<Mmap, std::io::Error> {
    let f = File::open(&fpath)?;
    let mmap = unsafe { Mmap::map(&f) }?;
    Ok(mmap)
}

/// UNSAFE
/// memory-map an existing BufReader
pub fn to_mmap(f: &File) -> Result<Mmap, std::io::Error> {
    let mmap = unsafe { Mmap::map(f) }?;
    Ok(mmap)
}

/// Convenience method to convert a file handle into a memory-mapped
/// fastq file reader
/// On Unix-alikes, will advise kernel of sequential memory access
pub fn mmapped_fq_reader<'a>(f: &'a File) -> Result<fastq::Reader<Cursor<Mmap>, fastq::Record>, std::io::Error> {
    let m = to_mmap(f)?;
    if cfg!(unix) {
        m.advise(Advice::Sequential)?;
    }
    let curs = Cursor::new(m);
    Ok(fastq::Reader::new(curs))
}

pub fn mmapped_cfq_reader<'a>(
    f: &'a File,
) -> Result<fastq::CheckedReader<Cursor<Mmap>, fastq::Record>, std::io::Error> {
    let m = to_mmap(f)?;
    if cfg!(unix) {
        m.advise(Advice::Sequential)?;
    }
    let curs = Cursor::new(m);
    Ok(fastq::CheckedReader::new(curs))
}

pub fn mmapped_fa_reader<'a>(f: &'a File) -> Result<fasta::Reader<Cursor<Mmap>, fasta::Record>, std::io::Error> {
    let m = to_mmap(f)?;
    if cfg!(unix) {
        m.advise(Advice::Sequential)?;
    }
    let curs = Cursor::new(m);
    Ok(fasta::Reader::new(curs))
}

pub fn mmapped_cfa_reader<'a>(f: &'a File) -> Result<fasta::CheckedReader<Cursor<Mmap>, fasta::Record>, std::io::Error> {
    let m = to_mmap(f)?;
    if cfg!(unix) {
        m.advise(Advice::Sequential)?;
    }
    let curs = Cursor::new(m);
    Ok(fasta::CheckedReader::new(curs))
}

pub fn is_gz(header: &[u8]) -> bool {
    (header[0] == GZB1) && (header[1] == GZB2)
}
