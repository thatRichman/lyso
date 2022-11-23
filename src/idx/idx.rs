use std::fs::File;
use std::io::prelude::*;

use crate::util;
use crate::io;

// NAME LENGTH OFFSET LINEBASES LINEWIDTH OFFSET


struct FileIndex<T> {
    index: T,
}

impl<F> FileIndex<F>
where F: io::fastq::Reader {
    
}

// Generic index trait, implement for each reader type
pub trait Index {
    type Kind;

    // load an existing index
    fn load_idx(&self, f: impl Read) -> FileIndex<T> {}
    
    // generate new index
    fn make_idx(&self) -> FileIndex<T> {}
    
    // write existing index to file
    fn write_idx(&self, idx: impl Write) -> File {}

    // get record of type Kind at index
    fn at<I>(&self, key: I) -> Option<Self::Kind> {}
 
}


#[cfg(test)]
mod tests {

    fn test_load_index(){}

    fn test_make_index(){}

    fn test_write_index(){}

    fn test_at(){}

}
