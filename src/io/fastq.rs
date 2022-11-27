use std::io;
use std::fmt;
use std::error::Error;
use std::path::Path;
use std::boxed::Box;
use std::io::Cursor;
use std::io::BufRead;
use std::iter::Iterator;
use std::iter::Map;
use itertools::Itertools;
use crate::io::mmap;
use crate::util::util;

#[derive(Debug)]
pub enum FastqError {
    ValidationError(&'static str),
    EofError,
    MissingId,
    TruncatedId,
    SeqQualMismatch,
    IoError(std::io::Error),
}

impl Error for FastqError {}

impl fmt::Display for FastqError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match self {
            Self::ValidationError(e) => write!(f, "Failed to validate record: {}", e),
            Self::EofError => write!(f, "Premature EOF encountered while parsing file"),
            Self::MissingId => write!(f, "Missing @ when attempting to read sequence ID"),
            Self::TruncatedId => write!(f, "Missing ID following @ symbol"),
            Self::SeqQualMismatch => write!(f, "Sequence length and quality length do not match"),
            Self::IoError(e) => write!(f, "IO Error: {}", e),
        }
    }
}

impl From<std::io::Error> for FastqError {
    fn from(e: std::io::Error) -> Self {
        Self::IoError(e)
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
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

   pub fn new() -> Self {
       Record {    
           id: "".to_string().into_boxed_str(),
           desc: "".to_string().into_boxed_str(),
           seq: "".to_string().into_boxed_str(),
           qual: "".to_string().into_boxed_str(),
       }
   }

	pub fn valid(&self) -> Result<(), FastqError> {
		

		if self.id.is_empty(){
			return Err(FastqError::ValidationError("id field cannot be empty"));
		}

		if !self.desc.is_ascii() {
			return Err(FastqError::ValidationError("Invalid character in description"));
		}
		
		if !self.seq.is_ascii() {
			return Err(FastqError::ValidationError("Invalid character in seq"));
		}

		if !self.qual.is_ascii() {
			return Err(FastqError::ValidationError("Invalid character in qual"));
		}
		
		if !self.seq.chars().all(|x| util::is_dna(x)){
			return Err(FastqError::ValidationError("invalid nucleotide in seq"));
		}
		
		if !self.seq.len() == self.qual.len() {
			return Err(FastqError::SeqQualMismatch);
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

    pub fn empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_empty() && self.seq.is_empty() && self.qual.is_empty()
    }

    pub fn clear(&mut self) {
        self.id = "".to_string().into_boxed_str();
        self.desc = "".to_string().into_boxed_str();
        self.seq = "".to_string().into_boxed_str();
        self.qual = "".to_string().into_boxed_str();
    }
 
}

pub struct Reader {
    inner: Cursor<Vec<u8>>,
    buffer: String,
}


impl Reader { 
	// new fastq reader from memory-mapped file
	pub fn from_mmapped_file<P: AsRef<Path> + std::fmt::Display>(fpath: P) -> Self {
        let mmap = mmap::to_mmap(&fpath);
		Reader {
            inner: Cursor::new(mmap),
            buffer: "".to_string(),
        }
        
	}
    
    pub fn read_record(&mut self, record: &mut Record) -> Result<(), FastqError>{ 
    
        // assume sequence is all on one line
        let mut vec_buf = Vec::with_capacity(1);

        record.clear();
        self.buffer.clear();
        
        match self.inner.read_line(&mut self.buffer) {
            Ok(0) if record.empty() => return Ok(()),  // EOF
            Ok(v) => self.buffer.trim_end().to_string(),
            Err(e) => {
                return Err(FastqError::IoError(e))
            },
        };

        if !self.buffer.starts_with("@") {
            return Err(FastqError::MissingId)
        }

        // assume all content after first whitespace is description
        let mut header = self.buffer[1..].trim_end().splitn(2, char::is_whitespace);
        match header.next() {
            Some(v) => { record.id = v.to_string().into_boxed_str()}
            None => return Err(FastqError::TruncatedId),
        }
        record.desc = header.next().unwrap_or("").to_string().into_boxed_str();
        self.buffer.clear();
        
        // read sequence content, track # of lines it takes up
        self.inner.read_line(&mut self.buffer)?;
        let mut seq_lines = 0;
        while !self.buffer.is_empty() && !self.buffer.starts_with("+") { 
            vec_buf.push(self.buffer.trim_end().to_string());
            self.buffer.clear();
            match self.inner.read_line(&mut self.buffer) {
                Ok(0) => return Err(FastqError::EofError),
                Ok(_) => (),
                Err(e) => return Err(FastqError::IoError(e)),
            }
            seq_lines += 1;
        }
        record.seq = vec_buf.join("").into_boxed_str(); 
        vec_buf.clear();
        
        // read the same # of lines for quality
        for _ in 0..seq_lines {
            self.buffer.clear();
            match self.inner.read_line(&mut self.buffer) {
                Ok(0) => {
                    return Err(FastqError::EofError)
                },
                Err(e) => {
                    return Err(FastqError::IoError(e))
                },
                Ok(v) => vec_buf.push(self.buffer.trim_end().to_string()),
            }
        }
        record.qual = vec_buf.join("").into_boxed_str();

        if record.qual.len() != record.seq.len() {
            return Err(FastqError::SeqQualMismatch)
        }

        Ok(())

    }
    
    pub fn position(&self) -> u64 {
        self.inner.position()
    }

    pub fn set_position(&mut self, pos: u64) {
        self.inner.set_position(pos);
    }

}

impl Iterator for Reader {
    type Item = Result<Record, FastqError>;
    
    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::new();
        match self.read_record(&mut record) {
            Ok(()) if record.empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(e) => Some(Err(e)),
        }
    }

}

#[cfg(test)]
mod tests {
    
    use crate::io::fastq::*;
    use std::path::PathBuf;
    
    fn init_path(s: &str) -> PathBuf {
        let mut test_data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        test_data_dir.push(s);
        test_data_dir
    }

    #[test]
    fn test_read_fq() {
        let fq_path = init_path("resources/test_data/test.fastq");
        let mut reader = Reader::from_mmapped_file(fq_path.as_path().display().to_string());
        assert!(reader.all(|x| x.expect("failed to read record").valid().is_ok()))
    }

    #[test]
    #[should_panic]
    fn test_bad_fq() {
        let fq_path = init_path("resources/test_data/trunc.fastq");
        let reader = Reader::from_mmapped_file(fq_path.as_path().display().to_string());
        for record in reader {
               match record {
                   Ok(_) => (), 
                   Err(e) => panic!("bad record: {}", e),
               } 
        }
    }
    

    #[test]
    fn test_get_fields() { 
        let fq_path = init_path("resources/test_data/test.fastq");
        let mut reader = Reader::from_mmapped_file(fq_path.as_path().display().to_string());
        let record = reader.next().unwrap().expect("bad record!");

        assert!(record.id == "R22092847.1.1".to_string().into_boxed_str());
        assert!(record.desc == "1 length=37".to_string().into_boxed_str());
        assert!(record.qual == "F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF".to_string().into_boxed_str());
        assert!(record.seq == "GNTTAAAGCACATAAAGACAAATCGCTCCAGGGCAAA".to_string().into_boxed_str()); 
    }

    #[test]
    fn test_fq_errors() {}
}
