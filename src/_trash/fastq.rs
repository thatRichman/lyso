use std::fs;
use std::error::Error;
use std::boxed::Box;
use crate::util;

pub const DEFAULT_CAPACITY: usize = 0;

// --- Begin FastQ Traits --- //


// --- End FastQ Traits --- //


// --- Begin FastQ Structs --- //

// Store fastq read information
// Use Box<str> because for most QC operations there is 
// no need for any of the fields to be resized, and
// we can perform masking for most mutating operations.
pub struct FQRead {
    id: Box<str>,
    read: Box<str>,
    qual: Box<str>,   
}

impl util::Validate for FQRead {
    fn seq_valid(&self) -> Result<bool, &'static str> {
        match self.read.chars().all(util::is_dna) {
            true => Ok(true),
            false => Err("Invalid character in read"),
        }
    }
    
    fn qual_valid(&self) -> Result<bool, &'static str> {
        return Ok(true);
    }
}

pub struct FQFile {
   file_path: String,
   reads: Vec<Box<str>>,
}

impl FQFile {
    pub fn new(file_path: String) -> Self {
        Self { 
            file_path: file_path,
            reads: Vec::with_capacity(DEFAULT_CAPACITY) 
        }
    }
}

// --- End FastQ Structs --- //


// --- Begin Tests --- //

#[cfg(test)]
mod test {
    use crate::fastq;
    use crate::util::Validate;
    use std::boxed::Box;

    #[test]
    fn test_fq_validate() {
        let id = String::from("test_id").into_boxed_str();
        let read = String::from("ATGCN").into_boxed_str();
        let qual = String::from("!#%$$^").into_boxed_str();
        let fq = fastq::FQRead{id: id,
                               read: read,
                               qual: qual
        }; 
        fq.valid().expect("valid FQRead is invalid?");
    }

    #[test]
    fn test_fq_invalid() { 
        let bad_read = String::from("ATGCNKLMOP").into_boxed_str();
        let bad_id = String::from("1").into_boxed_str();
        let bad_qual=String::from("$%&@!$@%").into_boxed_str();
        let fq_bad = fastq::FQRead{id: bad_id,
                                   read: bad_read,
                                   qual: bad_qual};
        match fq_bad.valid() {
            Ok(_val) => panic!("This should not be valid"),
            Err(err) => () 
        };
    }
}

// --- End Tests --- //

