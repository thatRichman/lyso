use std::fs::File;
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};

mod idx;
mod io;

mod util {
    pub mod util;
}

// mod idx {
//    pub mod idx;
// }

use crate::idx::FqIndexer;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Generate various file indices
    Faidx {
        /// Input file
        f_path: Option<PathBuf>,
    },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Faidx { f_path }) => {
            if let Some(p) = f_path.as_deref() {
                index_fastq(p);
            }
        }
        None => {}
    }

    fn index_fastq<P: AsRef<Path>>(fpath: P) {
        let mut in_file = File::open(&fpath).expect("unable to open file.");
        let mut buf_in = std::io::BufReader::new(&mut in_file);
        let fq_idxr = FqIndexer::new(&mut buf_in);
        let out_f = File::create(&"test.fai").unwrap();
        let mut buf_out = std::io::BufWriter::new(out_f);
        let mut idx_str;
        for idx in fq_idxr {
            idx_str = format!("{}\n", idx.unwrap());
            buf_out.write(idx_str.as_bytes()).unwrap();
        }
        buf_out.flush().unwrap();
    }
}
