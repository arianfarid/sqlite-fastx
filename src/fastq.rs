use crate::{
    SequenceCursor,
    fai::find_record_offset,
    filters::{CompareOp, ExecPlan, parse_plan},
    functions::{compute_gc, mean_quality, min_quality},
    reader::{ReadStrategy, SequenceReader, SequenceRecord},
};
use flate2::read::GzDecoder;
use seq_io::{fastq::*, policy::StdPolicy};
use sqlite3_ext::{Error, vtab::*, *};
use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
};

pub struct FastqSequenceReader {
    pub reader: Reader<Box<dyn Read>, StdPolicy>,
}
impl SequenceReader for FastqSequenceReader {
    type Record = OwnedRecord;

    fn next(&mut self) -> Option<Result<Self::Record>> {
        self.reader.next().map(|r| {
            r.map(|r| r.to_owned_record())
                .map_err(|e| sqlite3_ext::Error::from(e.to_string()))
        })
    }

    fn lookup_offset(fai_path: &str, id: &str) -> Option<u64> {
        use std::io::BufReader;
        let file = File::open(fai_path).ok()?;
        let mut reader = noodles_fastq::fai::io::Reader::new(BufReader::new(file));
        let mut buf = String::new();
        loop {
            buf.clear();
            match reader.read_record(&mut buf) {
                Ok(0) => return None,
                Ok(_) => {
                    if let Ok(record) = buf.parse::<noodles_fastq::fai::Record>() {
                        if record.name() == id {
                            return Some(record.sequence_offset());
                        }
                    }
                }
                Err(_) => return None,
            }
        }
    }
}

impl SequenceRecord for OwnedRecord {
    fn identifier_bytes(&self) -> &[u8] {
        seq_io::fastq::Record::id_bytes(self)
    }

    fn description_bytes(&self) -> Option<&[u8]> {
        seq_io::fastq::Record::desc_bytes(self)
    }

    fn sequence_bytes(&self) -> &[u8] {
        seq_io::fastq::Record::seq(self)
    }

    fn quality_bytes(&self) -> Option<&[u8]> {
        Some(seq_io::fastq::Record::qual(self))
    }
}

enum Columns {
    ID,
    Description,
    Sequence,
    Length,
    GCContent,
    Quality,
    MeanQuality,
    MinQuality,
    Filename,
}
impl TryFrom<i32> for Columns {
    type Error = ();

    fn try_from(value: i32) -> std::result::Result<Self, Self::Error> {
        match value {
            0 => Ok(Columns::ID),
            1 => Ok(Columns::Description),
            2 => Ok(Columns::Sequence),
            3 => Ok(Columns::Length),
            4 => Ok(Columns::GCContent),
            5 => Ok(Columns::Quality),
            6 => Ok(Columns::MeanQuality),
            7 => Ok(Columns::MinQuality),
            8 => Ok(Columns::Filename),
            _ => Err(()),
        }
    }
}

#[sqlite3_ext_vtab(StandardModule)]
pub struct FastqModule {
    filename: Option<String>,
    fai_path: Option<String>,
}
impl CreateVTab<'_> for FastqModule {
    fn create(
        _db: &'_ VTabConnection,
        _aux: &'_ Self::Aux,
        args: &[&str],
    ) -> Result<(String, Self)> {
        let filename = args
            .get(3)
            .map(|s| {
                let s = s.trim();
                let s = s.strip_prefix("filename=").unwrap_or(s);
                s.trim_matches('\'').to_string()
            })
            .unwrap();
        let schema = "CREATE TABLE x(
                id TEXT,
                description TEXT,
                sequence TEXT,
                length INTEGER,
                gc_content REAL,
                quality TEXT,
                mean_quality REAL,
                min_quality INTEGER,
                filename TEXT HIDDEN
            )";
        let fai_path = format!("{}.fai", filename);
        let fai_path = if std::path::Path::new(&fai_path).exists() {
            Some(fai_path)
        } else {
            None
        };
        Ok((
            schema.to_owned(),
            FastqModule {
                filename: Some(filename),
                fai_path,
            },
        ))
    }

    fn destroy(self) -> DisconnectResult<Self> {
        Ok(())
    }
}
impl VTab<'_> for FastqModule {
    type Aux = ();
    type Cursor = FastqCursor;
    fn connect(_db: &VTabConnection, _aux: &Self::Aux, args: &[&str]) -> Result<(String, Self)> {
        let filename = args
            .get(3)
            .map(|s| {
                let s = s.trim();
                let s = s.strip_prefix("filename=").unwrap_or(s);
                s.trim_matches('\'').to_string()
            })
            .unwrap();
        let schema = "CREATE TABLE x(
                id TEXT,
                description TEXT,
                sequence TEXT,
                length INTEGER,
                gc_content REAL,
                quality TEXT,
                mean_quality REAL,
                min_quality INTEGER,
                filename TEXT HIDDEN
            )";
        let fai_path = format!("{}.fai", filename);
        let fai_path = if std::path::Path::new(&fai_path).exists() {
            Some(fai_path)
        } else {
            None
        };
        Ok((
            schema.to_owned(),
            FastqModule {
                filename: Some(filename),
                fai_path,
            },
        ))
    }
    fn best_index(&self, index_info: &mut IndexInfo) -> Result<()> {
        let mut usable = vec![];
        for (i, constraint) in index_info.constraints().enumerate() {
            if constraint.usable() {
                match Columns::try_from(constraint.column())
                    .map_err(|_| Error::from("column index out of range"))?
                {
                    Columns::ID => match constraint.op() {
                        ConstraintOp::Like => usable.push((i, ("id", constraint.op()))),
                        ConstraintOp::Eq => usable.push((i, ("id", constraint.op()))),
                        _ => {}
                    },
                    #[allow(clippy::single_match)]
                    Columns::Description => match constraint.op() {
                        ConstraintOp::Like => usable.push((i, ("description", constraint.op()))),
                        _ => {}
                    },
                    #[allow(clippy::single_match)]
                    Columns::Sequence => match constraint.op() {
                        ConstraintOp::Like => usable.push((i, ("sequence", constraint.op()))),
                        _ => {}
                    },
                    Columns::Length => match constraint.op() {
                        ConstraintOp::GT
                        | ConstraintOp::GE
                        | ConstraintOp::LT
                        | ConstraintOp::LE
                        | ConstraintOp::Eq => {
                            usable.push((i, ("length", constraint.op())));
                        }
                        _ => {}
                    },
                    Columns::GCContent => match constraint.op() {
                        ConstraintOp::GT
                        | ConstraintOp::GE
                        | ConstraintOp::LT
                        | ConstraintOp::LE
                        | ConstraintOp::Eq => {
                            usable.push((i, ("gc_content", constraint.op())));
                        }
                        _ => {}
                    },
                    Columns::MeanQuality => match constraint.op() {
                        ConstraintOp::GT
                        | ConstraintOp::GE
                        | ConstraintOp::LT
                        | ConstraintOp::LE
                        | ConstraintOp::Eq => {
                            usable.push((i, ("mean_quality", constraint.op())));
                        }
                        _ => {}
                    },
                    Columns::MinQuality => match constraint.op() {
                        ConstraintOp::GT
                        | ConstraintOp::GE
                        | ConstraintOp::LT
                        | ConstraintOp::LE
                        | ConstraintOp::Eq => {
                            usable.push((i, ("min_quality", constraint.op())));
                        }
                        _ => {}
                    },
                    _ => {} // No op
                }
            }
        }
        if !usable.is_empty() {
            let mut constraints: Vec<_> = index_info.constraints().collect();

            let descriptor = usable
                .iter()
                .enumerate()
                .map(|(i, c)| {
                    constraints[c.0].set_argv_index(Some(i as u32));
                    constraints[c.0].set_omit(true);
                    let op_str = match c.1.1 {
                        ConstraintOp::GT => CompareOp::Gt.as_str(),
                        ConstraintOp::GE => CompareOp::Ge.as_str(),
                        ConstraintOp::LT => CompareOp::Lt.as_str(),
                        ConstraintOp::LE => CompareOp::Le.as_str(),
                        ConstraintOp::Eq => CompareOp::Eq.as_str(),
                        ConstraintOp::Like => "Like",
                        _ => "Scan",
                    };
                    let col_str = c.1.0;
                    [col_str, op_str].join(":")
                })
                .collect::<Vec<_>>()
                .join(",");
            index_info.set_index_str(Some(descriptor.as_str()))?;
        }
        index_info.set_estimated_cost(1000.0);
        Ok(())
    }
    fn open(&self) -> Result<Self::Cursor> {
        Ok(FastqCursor {
            plan: ExecPlan::new(),
            fallback_filename: self.filename.clone(),
            fai_path: self.fai_path.clone(),
            reader: None,
            current: None,
            rowid: 0,
            done: false,
            exit_early: false,
        })
    }
}
impl VTabCursor for SequenceCursor<FastqSequenceReader> {
    fn filter(
        &mut self,
        _index_num: i32,
        index_str: Option<&str>,
        args: &mut [&mut ValueRef],
    ) -> Result<()> {
        let strategy = self.determine_strategy(index_str, args)?;
        self.plan = parse_plan(index_str, args)?;

        let path = if let Some(ref f) = self.fallback_filename {
            f.clone()
        } else {
            return Err("filename constraint required".into());
        };

        let reader: Box<dyn Read> = match strategy {
            ReadStrategy::Stream => {
                if path.ends_with(".gz") {
                    let file = File::open(&path)
                        .map_err(|e| Error::from(format!("Cannot open '{}': {}", path, e)))?;
                    Box::new(GzDecoder::new(file))
                } else {
                    let file = File::open(&path)
                        .map_err(|e| Error::from(format!("Cannot open '{}': {}", path, e)))?;
                    Box::new(file)
                }
            }
            ReadStrategy::SeekToOffset(offset) => {
                let mut file = File::open(&path)
                    .map_err(|e| Error::from(format!("Cannot open '{}': {}", path, e)))?;

                let record_offset = find_record_offset(&mut file, offset).map_err(Error::from)?;

                file.seek(SeekFrom::Start(record_offset))
                    .map_err(|e| Error::from(format!("Seek failed: {}", e)))?;

                Box::new(file)
            }
        };
        self.reader = Some(FastqSequenceReader {
            reader: seq_io::fastq::Reader::new(reader),
        });
        self.rowid = 0;
        self.done = false;
        self.current = None;
        self.exit_early = false;
        self.next()
    }

    fn next(&mut self) -> Result<()> {
        let reader = self.reader.as_mut().ok_or("reader not initialized")?;
        loop {
            if self.exit_early {
                self.done = true;
                self.current = None;
                return Ok(());
            }
            match reader.next() {
                Some(Ok(record)) => {
                    if self.plan.eval(&record) {
                        self.current = Some(record);
                        self.rowid += 1;
                        if self.plan.unique {
                            self.exit_early = true;
                        }
                        return Ok(());
                    }
                }
                Some(Err(_)) => {
                    self.done = true;
                    return Ok(());
                }
                None => {
                    self.done = true;
                    self.current = None;
                    return Ok(());
                }
            }
        }
    }

    fn eof(&mut self) -> bool {
        self.done
    }

    fn column(&mut self, idx: usize, context: &ColumnContext) -> Result<()> {
        if let Some(record) = &self.current {
            match Columns::try_from(idx as i32)
                .map_err(|_| Error::from("column index out of range"))?
            {
                Columns::ID => context
                    .set_result(String::from_utf8_lossy(record.identifier_bytes()).to_string())?,
                Columns::Description => context.set_result(
                    record
                        .description_bytes()
                        .map(|d| String::from_utf8_lossy(d).to_string())
                        .unwrap_or_default(),
                )?,
                Columns::Sequence => context
                    .set_result(String::from_utf8_lossy(record.sequence_bytes()).to_string())?,
                Columns::Length => context.set_result(record.sequence_bytes().len() as i64)?,
                Columns::GCContent => context.set_result(compute_gc(record.sequence_bytes()))?,
                Columns::Quality => context.set_result(
                    record
                        .quality_bytes()
                        .map(|d| String::from_utf8_lossy(d).to_string())
                        .unwrap_or_default(),
                )?,
                Columns::MeanQuality => {
                    context.set_result(mean_quality(record.quality_bytes().unwrap_or_default()))?
                }
                Columns::MinQuality => {
                    context.set_result(min_quality(record.quality_bytes().unwrap_or_default()))?
                }
                Columns::Filename => {
                    context.set_result(self.fallback_filename.clone().unwrap_or_default())?
                }
            }
        }
        Ok(())
    }

    fn rowid(&mut self) -> Result<i64> {
        Ok(self.rowid)
    }
}

///Cursor for parsing FASTA files.
pub type FastqCursor = SequenceCursor<FastqSequenceReader>;
