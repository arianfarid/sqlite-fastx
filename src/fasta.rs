use crate::{
    SequenceCursor,
    fai::find_record_offset,
    filters::{CompareOp, ExecPlan, parse_plan},
    functions::compute_gc,
    reader::{ReadStrategy, SequenceReader, SequenceRecord},
};
use flate2::read::GzDecoder;
use seq_io::{fasta::*, policy::StdPolicy};
use sqlite3_ext::{Error, vtab::*, *};
use std::{
    fs::File,
    io::{BufReader, Read, Seek, SeekFrom},
};

pub struct FastaSequenceReader {
    pub reader: Reader<Box<dyn Read>, StdPolicy>,
}

impl SequenceReader for FastaSequenceReader {
    type Record = OwnedRecord;

    fn next(&mut self) -> Option<Result<Self::Record>> {
        self.reader.next().map(|r| {
            r.map(|r| r.to_owned_record())
                .map_err(|e| sqlite3_ext::Error::from(e.to_string()))
        })
    }

    fn lookup_offset(fai_path: &str, id: &str) -> Option<u64> {
        let index = match noodles_fasta::fai::fs::read(fai_path) {
            Ok(index) => index,
            Err(_) => return None,
        };
        let region = noodles_core::Region::new(id, ..);
        index.query(&region).ok()
    }
}
impl SequenceRecord for OwnedRecord {
    fn identifier_bytes(&self) -> &[u8] {
        seq_io::fasta::Record::id_bytes(self)
    }

    fn description_bytes(&self) -> Option<&[u8]> {
        seq_io::fasta::Record::desc_bytes(self)
    }

    fn sequence_bytes(&self) -> &[u8] {
        seq_io::fasta::Record::seq(self)
    }

    fn quality_bytes(&self) -> Option<&[u8]> {
        None
    }
}

enum Columns {
    ID,
    Description,
    Sequence,
    Length,
    GCContent,
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
            5 => Ok(Columns::Filename),
            _ => Err(()),
        }
    }
}

#[sqlite3_ext_vtab(StandardModule)]
pub struct FastaModule {
    filename: Option<String>,
    fai_path: Option<String>,
    is_bgzf: bool,
    indexed_columns: Vec<String>,
    index_fresh: bool,
}
const ALLOWED_RECORD_COLUMNS: &[&str] = &["id", "description", "length", "gc_content"];
const DEFAULT_RECORD_INDEXES: &[&str] = &["id", "length", "gc_content"];

impl CreateVTab<'_> for FastaModule {
    const SHADOW_NAMES: &'static [&'static str] = &["meta", "records"];
    fn create(
        db: &'_ VTabConnection,
        _aux: &'_ Self::Aux,
        args: &[&str],
    ) -> Result<(String, Self)> {
        let table_name = args.get(2).ok_or("missing table name")?;
        let filename = args
            .get(3)
            .map(|s| {
                let s = s.trim();
                let s = s.strip_prefix("filename=").unwrap_or(s);
                s.trim_matches('\'').to_string()
            })
            .unwrap();
        // Some kind of value, implicitly false
        let records_index = args
            .iter()
            .find(|arg| arg.starts_with("record_index="))
            .map_or("false", |arg| {
                arg.trim()
                    .strip_prefix("record_index=")
                    .unwrap_or(arg)
                    .trim_matches('\'')
            });
        let schema = "CREATE TABLE x(
                id TEXT,
                description TEXT,
                sequence TEXT,
                length INTEGER,
                gc_content REAL,
                filename TEXT HIDDEN
            )";

        let fai_path = format!("{}.fai", filename);
        let fai_path = if std::path::Path::new(&fai_path).exists() {
            Some(fai_path)
        } else {
            None
        };
        let is_bgzf = {
            if fai_path.is_none() {
                false
            } else {
                filename.ends_with(".gz")
            }
        };

        db.execute(
            &format!(
                "CREATE TABLE {table_name}_meta
        (source_file TEXT, indexed_columns TEXT, file_mtime INTEGER, file_size INTEGER, built_at INTEGER)"
            ),
            (),
        )?;
        db.execute(
            &format!(
                "CREATE TABLE {table_name}_records (id TEXT, description TEXT, length INTEGER, gc_content REAL, header_offset INTEGER)"
            ),
            (),
        )?;

        let metadata = std::fs::metadata(&filename).ok();
        let file_size = metadata.as_ref().map(|m| m.len() as i64);
        let file_mtime = metadata
            .as_ref()
            .and_then(|m| m.modified().ok())
            .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|d| d.as_secs() as i64);
        let built_at = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs() as i64)
            .unwrap_or(0);

        let indexed_columns: String = if records_index == "false" {
            "".to_string()
        } else if records_index == "true" {
            DEFAULT_RECORD_INDEXES.join(",")
        } else {
            let cols = records_index.replace(' ', "");
            if cols.split(',').all(|col| ALLOWED_RECORD_COLUMNS.contains(&col)) {
                cols
            } else {
                return Err(Error::from("Malformed".to_string()));
            }
        };

        let ico: Vec<String> = if indexed_columns.is_empty() {
            Vec::new()
        } else {
            indexed_columns.split(',').map(str::to_string).collect()
        };

        db.execute(
            &format!(
                "INSERT INTO {table_name}_meta (source_file, indexed_columns, file_mtime, file_size, built_at)
        VALUES (?,?,?,?,?)"
            ),
            params![filename.as_str(), &indexed_columns, file_mtime, file_size, built_at],
        )?;

        if !ico.is_empty() {
            let mut create_records_stmt = db.prepare(&format!(
                "INSERT INTO {table_name}_records (id, description, length, gc_content, header_offset) VALUES (?,?,?,?,?)"
            ))?;
            let is_compressed = filename.ends_with(".gz");
            let inner: Box<dyn Read> = if is_compressed {
                let file = File::open(&filename)
                    .map_err(|e| Error::from(format!("Cannot open '{}': {}", filename, e)))?;
                Box::new(GzDecoder::new(BufReader::new(file)))
            } else {
                let file = File::open(&filename)
                    .map_err(|e| Error::from(format!("Cannot open '{}': {}", filename, e)))?;
                Box::new(BufReader::new(file))
            };
            let mut fasta_reader = seq_io::fasta::Reader::new(inner);

            while let Some(result) = fasta_reader.next() {
                let record = result
                    .map_err(|e| Error::from(e.to_string()))?
                    .to_owned_record();
                let header_offset: Option<i64> = if is_compressed {
                    None
                } else {
                    fasta_reader.position().map(|p| p.byte() as i64)
                };
                let id = String::from_utf8_lossy(record.id_bytes());
                let desc = record
                    .desc_bytes()
                    .map(String::from_utf8_lossy)
                    .unwrap_or_default();
                let gc = compute_gc(record.sequence_bytes());
                create_records_stmt.execute(params!(
                    id.as_ref(),
                    desc.as_ref(),
                    record.seq().len() as i64,
                    gc,
                    header_offset
                ))?;
            }

            for col in &ico {
                db.execute(
                    &format!(
                        "CREATE INDEX {table_name}_records_{col}_idx ON {table_name}_records ({col});"
                    ),
                    (),
                )?;
            }
        }

        Ok((
            schema.to_owned(),
            FastaModule {
                filename: Some(filename),
                fai_path,
                is_bgzf,
                indexed_columns: ico,
                index_fresh: true,
            },
        ))
    }

    fn destroy(self) -> DisconnectResult<Self> {
        Ok(())
    }
}
impl VTab<'_> for FastaModule {
    type Aux = ();
    type Cursor = FastaCursor;
    fn connect(db: &VTabConnection, _aux: &Self::Aux, args: &[&str]) -> Result<(String, Self)> {
        let table_name = args.get(2).ok_or("missing table name")?;
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
                filename TEXT HIDDEN
            )";

        let fai_path = format!("{}.fai", filename);
        let fai_path = if std::path::Path::new(&fai_path).exists() {
            Some(fai_path)
        } else {
            None
        };
        let is_bgzf = fai_path.is_some() && filename.ends_with(".gz");

        let metadata = std::fs::metadata(&filename).ok();
        let file_size = metadata.as_ref().map(|m| m.len() as i64);
        let file_mtime = metadata
            .as_ref()
            .and_then(|m| m.modified().ok())
            .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|d| d.as_secs() as i64);
        let (index_column_row, built_mtime, built_size) = db.query_row(
            &format!("SELECT indexed_columns, file_mtime, file_size FROM {table_name}_meta"),
            (),
            |row| {
                let cols = row[0].get_str()?.to_string();
                let mtime = (!row[1].is_null()).then(|| row[1].get_i64());
                let size = (!row[2].is_null()).then(|| row[2].get_i64());
                Ok((cols, mtime, size))
            },
        )?;
        let indexed_columns: Vec<String> = if index_column_row.is_empty() {
            Vec::new()
        } else {
            index_column_row.split(',').map(str::to_string).collect()
        };
        let index_fresh = match (built_mtime, built_size, file_mtime, file_size) {
            (Some(bm), Some(bs), Some(m), Some(s)) => bm == m && bs == s,
            _ => false,
        };

        Ok((
            schema.to_owned(),
            FastaModule {
                filename: Some(filename),
                fai_path,
                is_bgzf,
                indexed_columns,
                index_fresh,
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
                    Columns::Sequence => {
                        #[allow(clippy::single_match)]
                        match constraint.op() {
                            ConstraintOp::Like => usable.push((i, ("sequence", constraint.op()))),
                            _ => {} //No op
                        }
                    }
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
        Ok(FastaCursor {
            plan: ExecPlan::new(),
            fallback_filename: self.filename.clone(),
            fai_path: self.fai_path.clone(),
            reader: None,
            current: None,
            rowid: 0,
            done: false,
            exit_early: false,
            is_bgzf: self.is_bgzf,
        })
    }
}

impl VTabCursor for SequenceCursor<FastaSequenceReader> {
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
                    if self.fai_path.is_some() {
                        Box::new(
                            noodles_bgzf::io::indexed_reader::Builder::default()
                                .build_from_path(&path)
                                .map_err(|e| {
                                    Error::from(format!("Cannot open '{}': {}", &path, e))
                                })?,
                        )
                    } else {
                        let file = File::open(&path)
                            .map_err(|e| Error::from(format!("Cannot open '{}': {}", &path, e)))?;
                        Box::new(GzDecoder::new(file))
                    }
                } else {
                    let file = File::open(&path)
                        .map_err(|e| Error::from(format!("Cannot open '{}': {}", path, e)))?;
                    Box::new(file)
                }
            }
            ReadStrategy::SeekToOffset(offset) => {
                let mut file = File::open(&path)
                    .map_err(|e| Error::from(format!("Cannot open '{}': {}", path, e)))?;
                if self.is_bgzf {
                    let mut indexed_reader = noodles_bgzf::io::indexed_reader::Builder::default()
                        .build_from_path(&path)
                        .map_err(|e| Error::from(format!("Cannot open '{}': {}", &path, e)))?;

                    let block_voffset = (offset >> 16) << 16;
                    indexed_reader
                        .seek(SeekFrom::Start(block_voffset))
                        .map_err(|e| Error::from(format!("Seek failed: {}", e)))?;
                    Box::new(indexed_reader)
                } else {
                    let record_offset =
                        find_record_offset(&mut file, offset).map_err(Error::from)?;
                    file.seek(SeekFrom::Start(record_offset))
                        .map_err(|e| Error::from(format!("Seek failed: {}", e)))?;

                    Box::new(file)
                }
            }
        };

        self.reader = Some(FastaSequenceReader {
            reader: seq_io::fasta::Reader::new(reader),
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
pub type FastaCursor = SequenceCursor<FastaSequenceReader>;
