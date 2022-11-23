use std::io;
use std::path::Path;
use std::boxed::Box;
use std::io::Cursor;
use std::io::BufRead;
use std::iter::Iterator;
use std::iter::Map;
use itertools::Itertools;

use crate::io::mmap;
use crate::util;

#[derive(Debug, Clone, Default)]
pub struct Record {
    id: Box<str>,
    desc: Box<str>,
    seq: Box<str>,
    qual: Box<str>,
}


// Records use Box<str> because they are
// not designed to need resizing, and the
// assumption is that there will be thousands
// to millions of them, so space saving is key
// Box will not actually allocate for empty str
impl Record {
   pub fn new(id: Box<str>, seq: Box<str>, desc:  Box<str>, qual: Box<str>) -> Self {
        
	    let r = Record {
			id: id, 
	    	desc: desc,
			seq: seq,
			qual: qual,
		};
        
        match r.valid() {
            Ok(_) => r,
            Err(e) => panic!("invalid record {:?}: {:?}", r, e),
        }
	}

	pub fn valid(&self) -> Result<(), &str> {
		

		if self.id.is_empty(){
			return Err("id field cannot be empty");
		}

		if !self.desc.is_ascii() {
			return Err("Invalid character in description");
		}
		
		if !self.seq.is_ascii() {
			return Err("Invalid character in seq");
		}

		if !self.qual.is_ascii() {
			return Err("Invalid character in qual");
		}
		
		if !self.seq.chars().all(|x| util::is_dna(x)){
			return Err("invalid nucleotide in seq");
		}
		
		if !self.seq.len() == self.qual.len() {
			return Err("seq and qual have differing lengths!");
		}
	
		Ok(())	
	}

	pub fn id(&self) -> &str {
		self.id.as_ref()
	}

	pub fn seq(&self) -> &str {
		self.seq.as_ref()
	}

	pub fn qual(&self) -> &str {
		self.qual.as_ref()
	}
 
}

pub struct Reader {
    inner: Cursor<Vec<u8>>,
}


impl Reader { 
	// new fastq reader from memory-mapped file
	pub fn from_mmapped_file<P: AsRef<Path> + std::fmt::Display>(fpath: P) -> Self {
        let mmap = mmap::to_mmap(&fpath);
		Reader {
            inner: Cursor::new(mmap)
        }
	}
    
    fn read_record(&mut self) -> Option<Vec<String>> {
        let mut lines: Vec<String> = Vec::with_capacity(4);

        for i in 0..4 {
            let mut buf = String::new();
            let bytes = self.inner.read_line(&mut buf);

            match bytes {
                Ok(0) if lines.len() == 0 => return None,  // EOF
                Ok(0) if lines.len() > 0 => {
                    panic!("premature EOF when reading record")
                },
                Err(e) => {
                    panic!("Encountered error reading record: {:?}", e)
                },
                _ => lines.push(buf
                        .trim()
                        .to_string()
                    ),
            };
        }

        Some(lines)
    }

}

impl Iterator for Reader {
    type Item = Record;
    
    fn next(&mut self) -> Option<Self::Item> {
        
        let mut rec = match self.read_record() {
            None => return None,
            Some(v) => v,
        };
        rec.reverse();
        Some(
            Record::new(
                rec
                    .pop()
                    .unwrap()
                    .into_boxed_str(),
                rec
                    .pop()
                    .unwrap()
                    .into_boxed_str(),
                rec
                    .pop()
                    .unwrap()
                    .into_boxed_str(),
                rec
                    .pop()
                    .unwrap()
                    .into_boxed_str(),
            )
        )
    }

}

#[cfg(test)]
#[should_panic]
fn test_valid() {
    let reader = Reader::from_mmapped_file("/home/spencer/projects/tetanus/test_data/test.fastq");
    let records = reader.collect::<Vec<Record>>();
    assert!(!records[0].valid().is_err(), "invalid fastq file!");
    
    let bad_read = "@test\nATGCX\n+S\nFFFFF".as_bytes();
    let bad_reader = Reader {
        inner: Cursor::new(bad_read.to_owned())
    };
    let bad_record = bad_reader.collect::<Vec<Record>>();
}

