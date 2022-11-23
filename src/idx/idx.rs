use std::fs::File;
use std::io::prelude::*;
use std::collections::HashMap;

use fxhash::FxHashMap;

use crate::util;
use crate::io;

struct FileIndex<'a, R: 'a, I> {
    reader: &'a R,
    index: I,
}

trait HasIndex {
    type Record;

    // load an existing index
    fn load_idx(&self, f: impl BufRead);
    
    // generate new index
    fn build_index(&self);
    
    // write existing index to file
    fn write_idx(&self) -> File;

    // get record at index
    fn at(&self, key: String) -> Self::Record;

}

impl<'a> HasIndex for FileIndex<'a, io::fastq::Reader, FxHashMap<String, [u32; 5]>> {
   
    // TODO this is naive and assumes a perfectly formed index
    fn load_idx(&self, f: impl BufRead) { 
        for entry in f.lines() {`
            let line = entry.unwrap().trim();
            let mut fields = line.split_whitespace();
            let key = fields.next().unwrap().to_string();
            let values = [
                fields.next().unwrap().parse::<u32>().unwrap(),
                fields.next().unwrap().parse::<u32>().unwrap(),
                fields.next().unwrap().parse::<u32>().unwrap(),
                fields.next().unwrap().parse::<u32>().unwrap(),
                fields.next().unwrap().parse::<u32>().unwrap(),
            ];

            self.index.insert(key, values); 
        }
    }

    fn build_idx(&mut self) {
        // reset cursor to beginning
        self.reader.set_position(0);
        assert_eq!(self.reader.position(), 0);
        

    } 
}


#[cfg(test)]
mod tests {

    fn test_load_index(){}

    fn test_make_index(){}

    fn test_wriyte_index(){}

    fn test_at(){}

}
