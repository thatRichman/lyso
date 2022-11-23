use std::io::Cursor;
use std::fs;
use std::path::Path;

use memmap2;

// UNSAFE
// given a file path, create a memory map to file contents
pub fn to_mmap<P: AsRef<Path> + std::fmt::Display> (fpath: P)-> Vec<u8>{
    let f = fs::File::open(&fpath).expect("Unable to open file");
    let mmap =  unsafe { 
        memmap2::Mmap::map(&f).expect(
            &format!("Error mapping file {}", fpath)
        )
    };
    mmap[..].to_vec()
}
