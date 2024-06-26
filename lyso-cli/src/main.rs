use std::fs::File;
use std::io::stdout;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::exit;

use bgzip;

use clap::{Parser, Subcommand};

use lyso_bam::reader::BamReader;
use lyso_bam::BamError;
use lyso_bam::Record as BamRecord;
use lyso_fasta::reader::FastaReader;
use lyso_fasta::FastaError;
use lyso_fastq::reader::FastqReader;
use lyso_fastq::FastqError;
use lyso_fastq::Record;

use std::time::Instant;

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
    View {
        f_path: Option<PathBuf>,
    },
    FaPrint {
        f_path: Option<PathBuf>,
    },
    FqPrint {
        f_path: Option<PathBuf>,
    },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Faidx { f_path }) => {
            if let Some(p) = f_path.as_deref() {
                unimplemented!();
            }
        }
        Some(Commands::View { f_path }) => {
            if let Some(p) = f_path.as_deref() {
                view_bam(p);
            }
        }
        Some(Commands::FaPrint { f_path }) => {
            if let Some(p) = f_path.as_deref() {
                test_read_fasta(p);
            }
        }
        Some(Commands::FqPrint { f_path }) => {
            if let Some(p) = f_path.as_deref() {
                test_read_fastq(p);
            }
        }
        None => {}
    }

    fn test_read_fasta<P: AsRef<Path>>(fpath: P) {
        let mut in_file = File::open(&fpath).expect("unable to open file.");
        let mut buf_in = std::io::BufReader::new(&mut in_file);
        let stdout = stdout();
        let mut handle = stdout.lock();
        let fa_reader = FastaReader::new(&mut buf_in);
        let now = Instant::now();
        let reads = fa_reader.collect::<Vec<Result<lyso_fasta::Record, FastaError>>>();
        eprintln!("Read {} records in {:?}", reads.len(), now.elapsed());
        // for rec in fa_reader {
        //     if let Err(e) = writeln!(handle, "{}", rec.unwrap()) {
        //         match e.kind() {
        //             std::io::ErrorKind::BrokenPipe => exit(141),
        //             _ => panic!("{e}"),
        //         }
        //     }
        // }
    }

    fn test_read_fastq<P: AsRef<Path>>(fpath: P) {
        let mut in_file = File::open(&fpath).expect("unable to open file.");
        let mut buf_in = std::io::BufReader::new(&mut in_file);
        let stdout = stdout();
        let mut handle = stdout.lock();
        let fa_reader = FastqReader::new(&mut buf_in);
        let now = Instant::now();
        let reads = fa_reader.collect::<Vec<Result<Record, FastqError>>>();
        eprintln!("Read {} records in {:?}", reads.len(), now.elapsed());
        // for rec in fa_reader {
        //     if let Err(e) = writeln!(handle, "{}", rec.unwrap()) {
        //         match e.kind() {
        //             std::io::ErrorKind::BrokenPipe => exit(141),
        //             _ => panic!("{e}"),
        //         }
        //     }
        // }
    }

    // fn index_fastq<P: AsRef<Path>>(fpath: P) {
    //     let mut in_file = File::open(&fpath).expect("unable to open file.");
    //     let mut buf_in = std::io::BufReader::new(&mut in_file);
    //     let fq_idxr = FastqIndexer::new(&mut buf_in);
    //     let out_f = File::create(&"test.fai").unwrap();
    //     let mut buf_out = std::io::BufWriter::new(out_f);
    //     let mut idx_str;
    //     for idx in fq_idxr {
    //         idx_str = format!("{}\n", idx.unwrap());
    //         buf_out.write(idx_str.as_bytes()).unwrap();
    //     }
    //     buf_out.flush().unwrap();
    // }

    fn view_bam<P: AsRef<Path>>(fpath: P) {
        let in_file = File::open(&fpath).expect("unable to open file.");
        let gunzip_in = bgzip::read::BGZFReader::new(in_file).unwrap();

        // automatically consume header and refs
        let bam_reader = BamReader::new(gunzip_in);
        let stdout = stdout();
        let mut handle = stdout.lock();
        //read alignments
        for rec in bam_reader.into_iter() {
            if let Err(e) = writeln!(handle, "{}", rec.unwrap()) {
                match e.kind() {
                    std::io::ErrorKind::BrokenPipe => exit(141),
                    _ => panic!("{e}"),
                }
            }
        }
    }
}
