use seq_io::{fasta::*, policy::StdPolicy};
use sqlite3_ext::{vtab::*, *};
use std::{fs::File, iter::Copied};

// ExecPlan is reserved for future optmizations
// Scan - full sequential FASTA scans
// IndexScan - planned
enum ExecPlan {
    Scan,
    IndexScan,
}

///Cursor for parsing FASTA files.
struct FastaCursor {
    plan: ExecPlan,
    fallback_filename: Option<String>,
    min_length: Option<i64>,
    max_length: Option<i64>,
    reader: Option<Reader<File, StdPolicy>>,
    current: Option<OwnedRecord>,
    rowid: i64,
    done: bool,
}
impl VTabCursor for FastaCursor {
    fn filter(
        &mut self,
        _index_num: i32,
        index_str: Option<&str>,
        args: &mut [&mut ValueRef],
    ) -> Result<()> {
        self.min_length = None;
        self.max_length = None;

        if !args.is_empty() {
            let val = args[0].get_i64();
            match index_str {
                Some("gt") => self.min_length = Some(val),
                Some("ge") => self.min_length = Some(val - 1),
                Some("lt") => self.max_length = Some(val),
                Some("le") => self.max_length = Some(val + 1),
                Some("eq") => {
                    self.min_length = Some(val - 1);
                    self.max_length = Some(val + 1);
                }
                _ => {}
            }
        }

        let path = if let Some(ref f) = self.fallback_filename {
            f.clone()
        } else {
            return Err("filename constraint required".into());
        };

        let file = File::open(&path).unwrap();

        self.reader = Some(seq_io::fasta::Reader::new(file));
        self.rowid = 0;
        self.done = false;
        self.current = None;
        self.next()
    }
    fn next(&mut self) -> Result<()> {
        let reader = self.reader.as_mut().unwrap();
        loop {
            match reader.next() {
                Some(Ok(record)) => {
                    let len = record.seq().len() as i64;
                    let passes = self.min_length.map_or(true, |m| len > m)
                        && self.max_length.map_or(true, |m| len < m);
                    if passes {
                        self.current = Some(record.to_owned_record());
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
                0 => context.set_result(String::from_utf8_lossy(record.id_bytes()).to_string())?,
                1 => context.set_result(
                    record
                        .desc_bytes()
                        .map(|d| String::from_utf8_lossy(d).to_string())
                        .unwrap_or_default(),
                )?,
                2 => context.set_result(String::from_utf8_lossy(&record.seq()).to_string())?,
                3 => context.set_result(record.seq().len() as i64)?,
                _ => {}
            }
        }
        Ok(())
    }
    fn rowid(&mut self) -> Result<i64> {
        Ok(self.rowid)
    }
}
#[sqlite3_ext_vtab(EponymousModule)]
struct FastaModule {
    filename: Option<String>,
}
impl VTab<'_> for FastaModule {
    type Aux = ();
    type Cursor = FastaCursor;
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
                filename TEXT HIDDEN
            )";
        Ok((
            schema.to_owned(),
            FastaModule {
                filename: Some(filename),
            },
        ))
    }
    fn best_index(&self, index_info: &mut IndexInfo) -> Result<()> {
        let mut arg_index: Option<usize> = None;

        for (i, constraint) in index_info.constraints().enumerate() {
            if constraint.usable() && constraint.column() == 3 {
                match constraint.op() {
                    ConstraintOp::GT
                    | ConstraintOp::GE
                    | ConstraintOp::LT
                    | ConstraintOp::LE
                    | ConstraintOp::Eq => {
                        arg_index = Some(i);
                    }
                    _ => {}
                }
            }
        }

        let mut constraints: Vec<_> = index_info.constraints().collect();
        if let Some(i) = arg_index {
            let descriptor = match index_info.constraints().nth(i).unwrap().op() {
                ConstraintOp::GT => "gt",
                ConstraintOp::GE => "ge",
                ConstraintOp::LT => "lt",
                ConstraintOp::LE => "le",
                ConstraintOp::Eq => "eq",
                _ => "scan",
            };
            constraints[i].set_argv_index(Some(0));
            constraints[i].set_omit(true);
            index_info.set_index_str(Some(descriptor))?;
        }

        index_info.set_estimated_cost(1000.0);
        Ok(())
    }
    fn open(&self) -> Result<Self::Cursor> {
        Ok(FastaCursor {
            plan: ExecPlan::Scan,
            fallback_filename: self.filename.clone(),
            reader: None,
            current: None,
            rowid: 0,
            done: false,
            min_length: None,
            max_length: None,
        })
    }
}
#[sqlite3_ext_main]
fn init(db: &Connection) -> Result<()> {
    db.create_module("fasta", FastaModule::module(), ())?;
    Ok(())
}
