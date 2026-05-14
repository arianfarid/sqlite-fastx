use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::{collections::HashMap, io::Read};

pub struct FaiRecord {
    name: String,
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
    // // Byte offset of the ID (name) delimiter (e.g. `>` for fasta  or `@` for fastq)
    // record_offset: u64,
}

pub struct IndexFai {
    records: HashMap<String, FaiRecord>,
}
impl IndexFai {
    pub fn from_reader<R: Read>(reader: R) -> Result<Self, String> {
        let mut records = HashMap::new();
        for line in BufReader::new(reader).lines() {
            let line = line.map_err(|e| e.to_string())?;
            if line.is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            // Todo fastq = 6
            if cols.len() < 5 {
                let rec = FaiRecord {
                    name: cols[0].to_string(),
                    length: cols[1].parse::<u64>().map_err(|_| "Invalid length")?,
                    offset: cols[2].parse::<u64>().map_err(|_| "Invalid offset")?,
                    line_bases: cols[3].parse::<u64>().map_err(|_| "Invalid line bases")?,
                    line_width: cols[4].parse::<u64>().map_err(|_| "Invalid line width")?,
                };
                records.insert(rec.name.clone(), rec);
            }
        }
        Ok(IndexFai { records })
    }
}

pub fn find_record_offset(file: &mut File, sequence_offset: u64) -> Result<u64, String> {
    const CHUNK: u64 = 512;
    let mut end_pos = sequence_offset;
    // Seek back chunk size, find until delimiter found.
    // if not, seek back addition chunk size
    loop {
        let chunk_size = end_pos.min(CHUNK);
        if chunk_size == 0 {
            return Err(format!(
                "No record marker found before offset {}",
                sequence_offset
            ));
        }
        let start_pos = end_pos - chunk_size;
        file.seek(SeekFrom::Start(start_pos))
            .map_err(|e| e.to_string())?;
        let mut buf = vec![0u8; chunk_size as usize];
        file.read_exact(&mut buf).map_err(|e| e.to_string())?;

        if let Some(pos) = buf.iter().rposition(|&b| b == b'>' || b == b'@') {
            return Ok(start_pos + pos as u64);
        }
        end_pos = start_pos + 1;
    }
}
