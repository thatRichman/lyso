#![allow(unused)]
use fxhash::{FxHashMap, FxHasher};
use std::collections::HashMap;
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{self, Seek, SeekFrom};
use std::io::{prelude::*, ErrorKind};
use std::marker::PhantomData;
use std::slice::IterMut;
use std::{default, fmt};

use crate::io::fastq::{self, FastqRecord};
use crate::io::fasta::{self, FastaRecord};
use crate::util;

// ****************************************** //
//               Fastq Indexing               //
// ****************************************** //

pub struct FastqIndex {
    inner: FxHashMap<String, FastqIndexEntry>,
}

impl FastqIndex {
    pub fn new() -> Self {
        FastqIndex {
            inner: FxHashMap::default(),
        }
    }

    pub fn from_entries<I>(entries: I) -> Self
    where
        I: Iterator<Item = FastqIndexEntry>,
    {
        let mut idx = Self::new();
        for e in entries {
            idx.inner.insert(e.name.clone(), e);
        }
        idx
    }

    pub fn from_fasta_file<F: BufRead + Seek>(fasta: &mut F) -> Self {
        let idxr = FastqIndexer::new(fasta);
        idxr.into()
    }

    pub fn read_index(&mut self, handle: &mut impl BufRead) -> Result<(), std::io::Error> {
        for line in handle.lines() {
            match line? {
                l => {
                    let fields = l.split('\t').collect::<Vec<&str>>();
                    if fields.len() != 6 {
                        return Err(std::io::Error::new(
                            ErrorKind::InvalidData,
                            "malformed index",
                        ));
                    }
                    self.inner.insert(
                        String::from(fields[0]),
                        FastqIndexEntry {
                            name: String::from(fields[0]),
                            offset: fields[2].parse::<u64>().unwrap(),
                            length: fields[1].parse::<u64>().unwrap(),
                            q_offset: fields[5].parse::<u64>().unwrap(),
                            linewidth: fields[4].parse::<u64>().unwrap(),
                            linebases: fields[3].parse::<u64>().unwrap(),
                        },
                    );
                }
            }
        }
        Ok(())
    }

    pub fn get(&self, id: &str) -> Option<&FastqIndexEntry> {
        self.inner.get(id)
    }

    pub fn inner(&self) -> &FxHashMap<String, FastqIndexEntry> {
        &self.inner
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct FastqIndexEntry {
    name: String,
    offset: u64,
    length: u64,
    q_offset: u64,
    linewidth: u64,
    linebases: u64,
}

impl fmt::Display for FastqIndexEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.name, self.length, self.offset, self.linebases, self.linewidth, self.q_offset
        )
    }
}

impl FastqIndexEntry {
    pub fn new() -> Self {
        FastqIndexEntry {
            name: "".into(),
            offset: 0,
            length: 0,
            q_offset: 0,
            linewidth: 0,
            linebases: 0,
        }
    }

    pub fn clear(&mut self) {
        self.name.clear();
        self.offset = 0;
        self.length = 0;
        self.q_offset = 0;
        self.linewidth = 0;
        self.linebases = 0;
    }

    pub fn empty(&self) -> bool {
        self.name.is_empty()
            && (*self.offset() == 0)
            && (*self.length() == 0)
            && (*self.q_offset() == 0)
            && (*self.linewidth() == 0)
            && (*self.linebases() == 0)
    }

    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn offset(&self) -> &u64 {
        &self.offset
    }

    pub fn length(&self) -> &u64 {
        &self.length
    }

    pub fn q_offset(&self) -> &u64 {
        &self.q_offset
    }

    pub fn linewidth(&self) -> &u64 {
        &self.linewidth
    }

    pub fn linebases(&self) -> &u64 {
        &self.linebases
    }
}

pub struct FastqIndexer<'a, R: 'a> {
    handle: &'a mut R,
    buffer: String,
}

impl<'a, F> FastqIndexer<'a, F>
where
    F: BufRead + Seek,
{
    pub fn new(f: &'a mut F) -> Self {
        FastqIndexer {
            handle: f,
            buffer: "".into(),
        }
    }

    pub fn make_index(&mut self, record: &mut FastqIndexEntry) -> Result<(), fastq::FastqError> {
        match self.handle.read_line(&mut self.buffer) {
            Ok(0) if record.empty() => return Ok(()), // EOF
            Ok(_) => self.buffer.retain(|c| c != '\n'),
            Err(e) => return Err(fastq::FastqError::IoError(e)),
        };

        if self.buffer.is_empty() && record.empty() {
            // EOF
            return Ok(());
        }

        if !self.buffer.starts_with("@") {
            return Err(fastq::FastqError::MissingId);
        }

        // assume all content after first whitespace is description
        let mut header = self.buffer[1..].trim_end().splitn(2, char::is_whitespace);
        match header.next() {
            Some(v) => record.name = v.to_string(),
            None => return Err(fastq::FastqError::TruncatedId),
        }
        record.offset = self.handle.stream_position()? as u64;

        // read first sequence line
        // don't count newline for nbases
        self.buffer.clear();
        record.linewidth = self.handle.read_line(&mut self.buffer)? as u64;
        record.linebases = self.buffer.trim_end().len() as u64;

        let mut seq_lines = 0;
        while !self.buffer.is_empty() && !self.buffer.starts_with("+") {
            record.length += self.buffer.trim_end().len() as u64;
            self.buffer.clear();
            match self.handle.read_line(&mut self.buffer) {
                Ok(0) => return Err(fastq::FastqError::EofError),
                Ok(_) => (),
                Err(e) => return Err(fastq::FastqError::IoError(e)),
            }
        }

        record.q_offset = (self.handle.stream_position()?) as u64;
        let skip_to = record.q_offset + record.length + 1;

        self.buffer.clear();
        // skip to start of next record
        self.handle.seek(SeekFrom::Start(skip_to));

        Ok(())
    }
}

impl<'a, F> Iterator for FastqIndexer<'a, F>
where
    F: BufRead + Seek,
{
    type Item = Result<FastqIndexEntry, fastq::FastqError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = FastqIndexEntry::new();
        match FastqIndexer::<'a, F>::make_index(self, &mut record) {
            Ok(()) if record.empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(e) => Some(Err(e)),
        }
    }
}

impl<'a, F> Into<FastqIndex> for FastqIndexer<'a, F>
where
    F: BufRead + Seek,
{
    fn into(self) -> FastqIndex {
        FastqIndex::from_entries(self.into_iter().map(|x| x.unwrap()))
    }
}

pub struct IndexedFastq<'a, F, R> {
    index: &'a FastqIndex,
    handle: F,
    _dtype: PhantomData<R>,
}

impl<'a, F, R> IndexedFastq<'a, F, R>
where
    F: BufRead + Seek,
    R: FastqRecord,
{
    pub fn new(handle: F, index: &'a FastqIndex) -> Self {
        IndexedFastq {
            index,
            handle,
            _dtype: PhantomData,
        }
    }

    pub fn get(&mut self, id: &str, rec: &mut R) -> Result<(), std::io::Error> {
        if let Some(idx) = self.index.get(id) {
            rec.clear();
            rec.set_id(idx.name.clone());

            self.handle.seek(SeekFrom::Start(idx.offset))?;
            let nlines = idx.length / idx.linebases;
            let mut buf: Vec<u8> = vec![0 as u8; (nlines * idx.linewidth) as usize];
            self.handle.read_exact(&mut buf)?;
            buf.retain(|c| *c != b'\n');
            rec.set_seq(String::from_utf8_lossy(&buf).into_owned());

            self.handle.seek(SeekFrom::Start(idx.q_offset))?;
            self.handle.read_exact(&mut buf)?;
            buf.retain(|c| *c != b'\n');
            rec.set_qual(String::from_utf8_lossy(&buf).into_owned());
            return Ok(());
        }
        Err(std::io::Error::new(ErrorKind::NotFound, "id not found"))
    }
}

// ****************************************** //
//               Fasta Indexing               //
// ****************************************** //

pub struct FastaIndex {
    inner: FxHashMap<String, FastaIndexEntry>,
}

impl FastaIndex {
    pub fn new() -> Self {
        FastaIndex {
            inner: FxHashMap::default(),
        }
    }

    pub fn from_entries<I>(entries: I) -> Self
    where
        I: Iterator<Item = FastaIndexEntry>,
    {
        let mut idx = Self::new();
        for e in entries {
            idx.inner.insert(e.name.clone(), e);
        }
        idx
    }

    pub fn from_fasta_file<F: BufRead + Seek>(fasta: &mut F) -> Self {
        let idxr = FastaIndexer::new(fasta);
        idxr.into()
    }

    pub fn read_index(&mut self, handle: &mut impl BufRead) -> Result<(), std::io::Error> {
        for line in handle.lines() {
            match line? {
                l => {
                    let fields = l.split('\t').collect::<Vec<&str>>();
                    if fields.len() != 5 {
                        return Err(std::io::Error::new(
                            ErrorKind::InvalidData,
                            "malformed index",
                        ));
                    }
                    self.inner.insert(
                        String::from(fields[0]),
                        FastaIndexEntry {
                            name: String::from(fields[0]),
                            offset: fields[2].parse::<u64>().unwrap(),
                            length: fields[1].parse::<u64>().unwrap(),
                            linewidth: fields[4].parse::<u64>().unwrap(),
                            linebases: fields[3].parse::<u64>().unwrap(),
                        },
                    );
                }
            }
        }
        Ok(())
    }

    pub fn get(&self, id: &str) -> Option<&FastaIndexEntry> {
        self.inner.get(id)
    }

    pub fn inner(&self) -> &FxHashMap<String, FastaIndexEntry> {
        &self.inner
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct FastaIndexEntry {
    name: String,
    offset: u64,
    length: u64,
    linewidth: u64,
    linebases: u64,
}

impl fmt::Display for FastaIndexEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.name, self.length, self.offset, self.linebases, self.linewidth
        )
    }
}

impl FastaIndexEntry {
    pub fn new() -> Self {
        FastaIndexEntry {
            name: "".into(),
            offset: 0,
            length: 0,
            linewidth: 0,
            linebases: 0,
        }
    }

    pub fn clear(&mut self) {
        self.name.clear();
        self.offset = 0;
        self.length = 0;
        self.linewidth = 0;
        self.linebases = 0;
    }

    pub fn empty(&self) -> bool {
        self.name.is_empty()
            && (*self.offset() == 0)
            && (*self.length() == 0)
            && (*self.linewidth() == 0)
            && (*self.linebases() == 0)
    }

    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn offset(&self) -> &u64 {
        &self.offset
    }

    pub fn length(&self) -> &u64 {
        &self.length
    }

    pub fn linewidth(&self) -> &u64 {
        &self.linewidth
    }

    pub fn linebases(&self) -> &u64 {
        &self.linebases
    }
}

pub struct FastaIndexer<'a, R: 'a> {
    handle: &'a mut R,
    buffer: String,
}

impl<'a, F> FastaIndexer<'a, F>
where
    F: BufRead + Seek,
{
    pub fn new(f: &'a mut F) -> Self {
        FastaIndexer {
            handle: f,
            buffer: "".into(),
        }
    }

    pub fn make_index(&mut self, record: &mut FastaIndexEntry) -> Result<(), fasta::FastaError> {
        match self.handle.read_line(&mut self.buffer) {
            Ok(0) if record.empty() => return Ok(()), // EOF
            Ok(_) => self.buffer.retain(|c| c != '\n'),
            Err(e) => return Err(fasta::FastaError::IoError(e)),
        };

        if self.buffer.is_empty() && record.empty() {
            // EOF
            return Ok(());
        }

        if !self.buffer.starts_with("@") {
            return Err(fasta::FastaError::MissingId);
        }

        // assume all content after first whitespace is description
        let mut header = self.buffer[1..].trim_end().splitn(2, char::is_whitespace);
        match header.next() {
            Some(v) => record.name = v.to_string(),
            None => return Err(fasta::FastaError::TruncatedId),
        }
        record.offset = self.handle.stream_position()? as u64;

        // read first sequence line
        // don't count newline for nbases
        self.buffer.clear();
        record.linewidth = self.handle.read_line(&mut self.buffer)? as u64;
        record.linebases = self.buffer.trim_end().len() as u64;

        let mut seq_lines = 0;
        while !self.buffer.is_empty() && !self.buffer.starts_with("+") {
            record.length += self.buffer.trim_end().len() as u64;
            self.buffer.clear();
            match self.handle.read_line(&mut self.buffer) {
                Ok(0) => return Err(fasta::FastaError::EofError),
                Ok(_) => (),
                Err(e) => return Err(fasta::FastaError::IoError(e)),
            }
        }

        let skip_to = record.offset + record.length + 1;
        self.buffer.clear();
        // skip to start of next record
        self.handle.seek(SeekFrom::Start(skip_to));

        Ok(())
    }
}

impl<'a, F> Iterator for FastaIndexer<'a, F>
where
    F: BufRead + Seek,
{
    type Item = Result<FastaIndexEntry, fasta::FastaError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = FastaIndexEntry::new();
        match FastaIndexer::<'a, F>::make_index(self, &mut record) {
            Ok(()) if record.empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(e) => Some(Err(e)),
        }
    }
}

impl<'a, F> Into<FastaIndex> for FastaIndexer<'a, F>
where
    F: BufRead + Seek,
{
    fn into(self) -> FastaIndex {
        FastaIndex::from_entries(self.into_iter().map(|x| x.unwrap()))
    }
}

pub struct IndexedFasta<'a, F, R> {
    index: &'a FastaIndex,
    handle: F,
    _dtype: PhantomData<R>,
}

impl<'a, F, R> IndexedFasta<'a, F, R>
where
    F: BufRead + Seek,
    R: FastaRecord,
{
    pub fn new(handle: F, index: &'a FastaIndex) -> Self {
        IndexedFasta {
            index,
            handle,
            _dtype: PhantomData,
        }
    }

    pub fn get(&mut self, id: &str, rec: &mut R) -> Result<(), std::io::Error> {
        if let Some(idx) = self.index.get(id) {
            rec.clear();
            rec.set_id(idx.name.clone());

            self.handle.seek(SeekFrom::Start(idx.offset))?;
            let nlines = idx.length / idx.linebases;
            let mut buf: Vec<u8> = vec![0 as u8; (nlines * idx.linewidth) as usize];
            self.handle.read_exact(&mut buf)?;
            buf.retain(|c| *c != b'\n');
            rec.set_seq(String::from_utf8_lossy(&buf).into_owned());

            return Ok(());
        }
        Err(std::io::Error::new(ErrorKind::NotFound, "id not found"))
    }
}

#[cfg(test)]
mod tests {
    use crate::idx;
    use crate::io;
    use crate::util;

    fn test_make_index() {}
}
