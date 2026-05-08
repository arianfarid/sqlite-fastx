use crate::{
    SequenceCursor,
    filters::{CompareOp, ExecPlan, parse_plan},
    functions::compute_gc,
    reader::{SequenceReader, SequenceRecord},
};
use flate2::read::GzDecoder;
use seq_io::{fastq::*, policy::StdPolicy};
use sqlite3_ext::{Error, vtab::*, *};
use std::{fs::File, io::Read};

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

#[sqlite3_ext_vtab(EponymousModule)]
pub struct FastqModule {
    filename: Option<String>,
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
                filename TEXT HIDDEN
            )";
        Ok((
            schema.to_owned(),
            FastqModule {
                filename: Some(filename),
            },
        ))
    }
    fn best_index(&self, index_info: &mut IndexInfo) -> Result<()> {
        let mut usable = vec![];
        for (i, constraint) in index_info.constraints().enumerate() {
            if constraint.usable() {
                match constraint.column() {
                    // Sequence
                    2 => {
                        match constraint.op() {
                            ConstraintOp::Like => usable.push((i, ("sequence", constraint.op()))),
                            _ => {} //No op
                        }
                    }
                    // Length
                    3 => match constraint.op() {
                        ConstraintOp::GT
                        | ConstraintOp::GE
                        | ConstraintOp::LT
                        | ConstraintOp::LE
                        | ConstraintOp::Eq => {
                            usable.push((i, ("length", constraint.op())));
                        }
                        _ => {}
                    },
                    // GC Content
                    4 => match constraint.op() {
                        ConstraintOp::GT
                        | ConstraintOp::GE
                        | ConstraintOp::LT
                        | ConstraintOp::LE
                        | ConstraintOp::Eq => {
                            usable.push((i, ("gc_content", constraint.op())));
                        }
                        _ => {}
                    },
                    // Quality
                    // 5 => match
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
            reader: None,
            current: None,
            rowid: 0,
            done: false,
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
        self.plan = parse_plan(index_str, args)?;

        let path = if let Some(ref f) = self.fallback_filename {
            f.clone()
        } else {
            return Err("filename constraint required".into());
        };

        let reader: Box<dyn Read> = if path.ends_with(".gz") {
            let file = File::open(&path)
                .map_err(|e| Error::from(format!("Cannot open file '{}': {}", path, e)))?;
            Box::new(GzDecoder::new(file))
        } else {
            let file = File::open(&path)
                .map_err(|e| return Error::from(format!("Cannot open file '{}': {}", path, e)))?;
            Box::new(file)
        };
        self.reader = Some(FastqSequenceReader {
            reader: seq_io::fastq::Reader::new(reader),
        });
        self.rowid = 0;
        self.done = false;
        self.current = None;
        self.next()
    }

    fn next(&mut self) -> Result<()> {
        let reader = self
            .reader
            .as_mut()
            .ok_or_else(|| "reader not initialized")?;
        loop {
            match reader.next() {
                Some(Ok(record)) => {
                    if self.plan.matches(&record) {
                        self.current = Some(record.clone());
                        self.rowid += 1;
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
            match idx {
                0 => context
                    .set_result(String::from_utf8_lossy(record.identifier_bytes()).to_string())?,
                1 => context.set_result(
                    record
                        .description_bytes()
                        .map(|d| String::from_utf8_lossy(d).to_string())
                        .unwrap_or_default(),
                )?,
                2 => context
                    .set_result(String::from_utf8_lossy(&record.sequence_bytes()).to_string())?,
                3 => context.set_result(record.sequence_bytes().len() as i64)?,
                4 => context.set_result(compute_gc(record.sequence_bytes()))?,
                5 => context.set_result(
                    record
                        .quality_bytes()
                        .map(|d| String::from_utf8_lossy(d).to_string())
                        .unwrap_or_default(),
                )?,
                _ => {}
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
