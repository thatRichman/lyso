mod io {
    pub mod fastq;
    pub mod mmap;
}

mod util {
    pub mod util;
}

mod idx {
    pub mod idx;
}

use crate::io::fastq;

fn main() {
    
    let reader = fastq::Reader::from_mmapped_file("/home/spencer/projects/tetanus/test_data/test.fastq ");

}
