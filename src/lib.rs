use seq_io::{fasta::*, policy::StdPolicy};
use sqlite3_ext::{Error, function::FunctionOptions, vtab::*, *};
use std::fs::File;

enum Predicate {
    Length(LengthFilter),
    SequenceLike(SequenceFilter),
    // TODO: GC, Substring
}
impl Predicate {
    fn matches(&self, record: &RefRecord) -> bool {
        match self {
            Predicate::Length(f) => f.matches(record.seq().len() as i64),
            Predicate::SequenceLike(s) => s.like(record.seq()),
        }
    }
}
#[repr(i32)]
enum LengthOp {
    Gt,
    Ge,
    Lt,
    Le,
    Eq,
}
impl LengthOp {
    fn as_str(&self) -> &'static str {
        match self {
            LengthOp::Gt => "Gt",
            LengthOp::Ge => "Ge",
            LengthOp::Lt => "Lt",
            LengthOp::Le => "Le",
            LengthOp::Eq => "Eq",
        }
    }
}

enum SequenceOp {
    Contains,
    StartsWith,
    EndsWith,
    Eq,
}
impl SequenceOp {
    fn as_str(&self) -> &'static str {
        match self {
            SequenceOp::Contains => "Contains",
            SequenceOp::StartsWith => "StartsWith",
            SequenceOp::EndsWith => "EndsWith",
            SequenceOp::Eq => "Eq",
        }
    }
}
struct LengthFilter {
    op: LengthOp,
    value: i64,
}
impl LengthFilter {
    fn matches(&self, len: i64) -> bool {
        match self.op {
            LengthOp::Gt => len > self.value,
            LengthOp::Ge => len >= self.value,
            LengthOp::Lt => len < self.value,
            LengthOp::Le => len <= self.value,
            LengthOp::Eq => len == self.value,
        }
    }
}
struct SequenceFilter {
    op: SequenceOp,
    pattern: String,
}
impl SequenceFilter {
    fn like(&self, seq: &[u8]) -> bool {
        match self.op {
            SequenceOp::Contains => seq
                .windows(self.pattern.len())
                .any(|w| w == self.pattern.as_bytes()),
            SequenceOp::StartsWith => seq.starts_with(self.pattern.as_bytes()),
            SequenceOp::EndsWith => seq.ends_with(self.pattern.as_bytes()),
            SequenceOp::Eq => seq == self.pattern.as_bytes(),
        }
    }
}

struct ExecPlan {
    predicates: Vec<Predicate>,
}
impl ExecPlan {
    fn new() -> ExecPlan {
        ExecPlan { predicates: vec![] }
    }
    fn matches(&self, record: &RefRecord) -> bool {
        for pred in &self.predicates {
            if !pred.matches(record) {
                return false;
            }
        }
        true
    }
}

///Cursor for parsing FASTA files.
struct FastaCursor {
    plan: ExecPlan,
    fallback_filename: Option<String>,
    reader: Option<Reader<File, StdPolicy>>,
    current: Option<OwnedRecord>,
    rowid: i64,
    done: bool,
}
impl FastaCursor {
    fn parse_plan(index_str: Option<&str>, args: &mut [&mut ValueRef]) -> Result<ExecPlan> {
        let Some(descriptor) = index_str else {
            return Ok(ExecPlan { predicates: vec![] });
        };

        let mut predicates = vec![];
        for (i, op_str) in descriptor.split(',').enumerate() {
            if i >= args.len() {
                break;
            }
            let arg = &mut args[i];
            match op_str {
                "Gt" => predicates.push(Predicate::Length(LengthFilter {
                    op: LengthOp::Gt,
                    value: arg.get_i64(),
                })),
                "Ge" => predicates.push(Predicate::Length(LengthFilter {
                    op: LengthOp::Ge,
                    value: arg.get_i64(),
                })),
                "Lt" => predicates.push(Predicate::Length(LengthFilter {
                    op: LengthOp::Lt,
                    value: arg.get_i64(),
                })),
                "Le" => predicates.push(Predicate::Length(LengthFilter {
                    op: LengthOp::Le,
                    value: arg.get_i64(),
                })),
                "Eq" => predicates.push(Predicate::Length(LengthFilter {
                    op: LengthOp::Eq,
                    value: arg.get_i64(),
                })),
                "Like" => {
                    let raw = arg.get_str()?.to_string();
                    let starts_with_wild = raw.starts_with('%');
                    let ends_with_wild = raw.ends_with('%');
                    let pattern = raw.trim_matches('%').to_ascii_uppercase();
                    let op = match (starts_with_wild, ends_with_wild) {
                        (true, true) => SequenceOp::Contains,
                        (true, false) => SequenceOp::StartsWith,
                        (false, true) => SequenceOp::EndsWith,
                        (false, false) => SequenceOp::Eq,
                    };
                    predicates.push(Predicate::SequenceLike(SequenceFilter { op, pattern }))
                }
                _ => continue,
            };
        }
        Ok(ExecPlan { predicates })
    }
}

impl VTabCursor for FastaCursor {
    fn filter(
        &mut self,
        _index_num: i32,
        index_str: Option<&str>,
        args: &mut [&mut ValueRef],
    ) -> Result<()> {
        self.plan = Self::parse_plan(index_str, args)?;

        let path = if let Some(ref f) = self.fallback_filename {
            f.clone()
        } else {
            return Err("filename constraint required".into());
        };

        let file = File::open(&path)
            .map_err(|e| return Error::from(format!("Cannot open file '{}': {}", path, e)))?;

        self.reader = Some(seq_io::fasta::Reader::new(file));
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
        let mut usable = vec![];
        for (i, constraint) in index_info.constraints().enumerate() {
            // Length column constraints
            if constraint.usable() {
                match constraint.column() {
                    2 => {
                        match constraint.op() {
                            ConstraintOp::Like => usable.push((i, constraint.op())),
                            _ => {} //No op
                        }
                    }
                    3 => match constraint.op() {
                        ConstraintOp::GT
                        | ConstraintOp::GE
                        | ConstraintOp::LT
                        | ConstraintOp::LE
                        | ConstraintOp::Eq => {
                            usable.push((i, constraint.op()));
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
                    match c.1 {
                        ConstraintOp::GT => LengthOp::Gt.as_str(),
                        ConstraintOp::GE => LengthOp::Ge.as_str(),
                        ConstraintOp::LT => LengthOp::Lt.as_str(),
                        ConstraintOp::LE => LengthOp::Le.as_str(),
                        ConstraintOp::Eq => LengthOp::Eq.as_str(),
                        ConstraintOp::Like => "Like",
                        _ => "Scan",
                    }
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
            reader: None,
            current: None,
            rowid: 0,
            done: false,
        })
    }
}

fn compute_gc(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
        .count();
    gc as f64 / seq.len() as f64
}

#[sqlite3_ext_main]
pub fn init(db: &Connection) -> Result<()> {
    db.create_module("fasta", FastaModule::module(), ())?;
    db.create_scalar_function("gc_content", &FunctionOptions::default(), |ctx, args| {
        let seq = args[0].get_str()?;
        let gc = compute_gc(seq.as_bytes());
        ctx.set_result(gc)
    })?;
    Ok(())
}
