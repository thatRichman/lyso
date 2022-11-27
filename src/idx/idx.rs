use std::fs::File;
use std::io::prelude::*;
use std::collections::HashMap;

use fxhash::FxHashMap;

use crate::util;
use crate::io;

type FqIdx = [u64; 5];

struct FileIndex<'a, R: 'a, I> {
    reader: &'a mut R,
    index: I,
}

trait HasIndex {
    type Record;

    // load an existing index
    fn load_idx(&mut self, f: impl BufRead);
    
    // generate new index
    fn build_idx(&mut self);
    
    // write existing index to file
    fn write_idx(&self, f: impl Write);

    // get record at index
    fn at(&mut self, key: String);

}

impl<'a> HasIndex for FileIndex<'a, io::fastq::Reader, FxHashMap<String, FqIdx>> {
    
    type Record = io::fastq::Record;

    // TODO this is naive and assumes a perfectly formed index
    fn load_idx(&mut self, f: impl BufRead) { 
        for entry in f.lines() {
            let line = entry.unwrap();
            let line = line.trim();
            let mut fields = line.split_whitespace();
            let key = fields.next().unwrap().to_string();
            let values = [
                fields.next().unwrap().parse::<u64>().unwrap(),
                fields.next().unwrap().parse::<u64>().unwrap(),
                fields.next().unwrap().parse::<u64>().unwrap(),
                fields.next().unwrap().parse::<u64>().unwrap(),
                fields.next().unwrap().parse::<u64>().unwrap(),
            ];

            self.index.insert(key, values); 
        }
    }

    fn build_idx(&mut self) {
        // reset cursor to beginning
        self.reader.set_position(0);
    }

    fn write_idx(&self, f: impl Write) {
    
    }

    fn at(&mut self, key: String) {
        // get record offset
        // set reader position to offset
        // read record
    }
}


#[cfg(test)]
mod tests {
    use crate::io;
    use crate::util;
    use crate::idx;

    fn test_load_index(){
        
    }

    fn test_make_index(){}

    fn test_write_index(){}

    fn test_at(){}

}
