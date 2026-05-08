use seq_io::{fasta::*, policy::StdPolicy};
use sqlite3_ext::{Error, function::FunctionOptions, vtab::*, *};
use std::fs::File;

enum Predicate {
    Length(LengthFilter),
    SequenceLike(SequenceFilter),
    // TODO: GC, Substring
}
impl Predicate {
    fn matches<S: SequenceRecord>(&self, record: &S) -> bool {
        match self {
            Predicate::Length(f) => f.matches(record.sequence_bytes().len() as i64),
            Predicate::SequenceLike(s) => s.like(record.sequence_bytes()),
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
    fn matches<S: SequenceRecord>(&self, record: &S) -> bool {
        for pred in &self.predicates {
            if !pred.matches(record) {
                return false;
            }
        }
        true
    }
}

trait SequenceRecord: Clone {
    fn identifier_bytes(&self) -> &[u8];
    fn description_bytes(&self) -> Option<&[u8]>;
    fn sequence_bytes(&self) -> &[u8];
}
trait SequenceReader {
    type Record: SequenceRecord;
    fn next(&mut self) -> Option<Result<Self::Record>>;
}

struct FastaSequenceReader {
    reader: Reader<File, StdPolicy>,
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
}
impl SequenceReader for FastaSequenceReader {
    type Record = OwnedRecord;

    fn next(&mut self) -> Option<Result<Self::Record>> {
        self.reader.next().map(|r| {
            r.map(|r| r.to_owned_record())
                .map_err(|e| sqlite3_ext::Error::from(e.to_string()))
        })
    }
}

struct SequenceCursor<R: SequenceReader> {
    plan: ExecPlan,
    fallback_filename: Option<String>,
    reader: Option<R>,
    current: Option<R::Record>,
    rowid: i64,
    done: bool,
}
impl VTabCursor for SequenceCursor<FastaSequenceReader> {
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

        let file = File::open(&path)
            .map_err(|e| return Error::from(format!("Cannot open file '{}': {}", path, e)))?;

        self.reader = Some(FastaSequenceReader {
            reader: seq_io::fasta::Reader::new(file),
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

type FastaCursor = SequenceCursor<FastaSequenceReader>;
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

fn dna_to_rna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'T' => b'U',
            b't' => b'u',
            _ => b,
        })
        .collect()
}

fn rna_to_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'U' => b'T',
            b'u' => b't',
            _ => b,
        })
        .collect()
}

fn n_count(seq: &[u8]) -> i64 {
    seq.iter().filter(|&b| matches!(b, b'n' | b'N')).count() as i64
}

fn base_count(seq: &[u8], base: u8) -> sqlite3_ext::Result<i64> {
    if !matches!(
        base.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T' | b'U' | b'N'
    ) {
        return Err("Not a valid base".into());
    }
    let base = base.to_ascii_uppercase();
    Ok(seq
        .iter()
        .filter(|b| b.to_ascii_uppercase() == base)
        .count() as i64)
}

fn is_valid_dna(seq: &[u8]) -> bool {
    seq.iter()
        .all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}
fn is_valid_rna(seq: &[u8]) -> bool {
    seq.iter()
        .all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'U' | b'N'))
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'a' => b't',
            b'G' => b'C',
            b'g' => b'c',
            b'T' => b'A',
            b't' => b'a',
            b'C' => b'G',
            b'c' => b'g',
            b'U' => b'A',
            b'u' => b'a',
            _ => b,
        })
        .collect()
}

#[sqlite3_ext_main]
pub fn init(db: &Connection) -> Result<()> {
    db.create_module("fasta", FastaModule::module(), ())?;
    db.create_scalar_function(
        "gc_content",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let gc = compute_gc(seq.as_bytes());
            ctx.set_result(gc)
        },
    )?;
    db.create_scalar_function(
        "n_count",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let count = n_count(seq.as_bytes());
            ctx.set_result(count)
        },
    )?;
    db.create_scalar_function(
        "base_count",
        &FunctionOptions::default().set_n_args(2),
        |ctx, args| {
            let seq = args[0].get_str()?.to_string();
            let base_str = &args[1].get_str()?.to_string();
            let base = base_str
                .as_bytes()
                .first()
                .ok_or_else(|| "base_count requires a non-empty base argument")?;
            let count = base_count(seq.as_bytes(), *base)?;
            ctx.set_result(count)
        },
    )?;
    db.create_scalar_function(
        "to_rna",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let seq = dna_to_rna(seq.as_bytes());
            ctx.set_result(String::from_utf8_lossy(&seq).into_owned())
        },
    )?;
    db.create_scalar_function(
        "to_dna",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let seq = rna_to_dna(seq.as_bytes());
            ctx.set_result(String::from_utf8_lossy(&seq).into_owned())
        },
    )?;
    db.create_scalar_function(
        "reverse_complement",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let seq = reverse_complement(seq.as_bytes());
            ctx.set_result(String::from_utf8_lossy(&seq).into_owned())
        },
    )?;
    db.create_scalar_function(
        "is_valid_dna",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let valid = is_valid_dna(seq.as_bytes());
            ctx.set_result(valid)
        },
    )?;
    db.create_scalar_function(
        "is_valid_rna",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let valid = is_valid_rna(seq.as_bytes());
            ctx.set_result(valid)
        },
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    // compute_gc tests
    #[test]
    fn gc_empty() {
        assert_eq!(compute_gc(b""), 0.0);
    }

    #[test]
    fn gc_all_gc() {
        assert_eq!(compute_gc(b"GGCC"), 1.0);
        assert_eq!(compute_gc(b"ggcc"), 1.0);
    }

    #[test]
    fn gc_no_gc() {
        assert_eq!(compute_gc(b"AATT"), 0.0);
        assert_eq!(compute_gc(b"aatt"), 0.0);
    }

    #[test]
    fn gc() {
        assert_eq!(compute_gc(b"ACGT"), 0.5);
    }

    #[test]
    fn gc_with_ambiguous_bases() {
        assert_eq!(compute_gc(b"ATGCN"), 0.4);
    }

    // n_count
    #[test]
    fn ncount_empty() {
        assert_eq!(n_count(b""), 0);
    }

    #[test]
    fn ncount_all_n() {
        assert_eq!(n_count(b"NNNN"), 4);
    }

    #[test]
    fn ncount_no_n() {
        assert_eq!(n_count(b"ACGT"), 0);
    }

    #[test]
    fn ncount_3n() {
        assert_eq!(n_count(b"ACNGNTN"), 3);
    }

    #[test]
    fn ncount_lower_n() {
        assert_eq!(n_count(b"ACnGnTn"), 3);
    }

    // base_count
    #[test]
    fn base_count_basic() {
        assert_eq!(base_count(b"ACGT", b'A'), Ok(1));
        assert_eq!(base_count(b"AAAA", b'A'), Ok(4));
    }

    #[test]
    fn base_count_case_insensitive() {
        assert_eq!(base_count(b"AaCcGgTt", b'a'), Ok(2));
        assert_eq!(base_count(b"AaCcGgTt", b'A'), Ok(2));
    }

    #[test]
    fn base_count_invalid_base() {
        assert!(base_count(b"ACGT", b'Z').is_err());
    }

    #[test]
    fn base_count_empty() {
        assert_eq!(base_count(b"", b'A'), Ok(0));
    }

    //to_rna
    #[test]
    fn rna_converts_t_to_u() {
        assert_eq!(dna_to_rna(b"ACGT"), b"ACGU");
    }

    #[test]
    fn rna_lowercase_t() {
        assert_eq!(dna_to_rna(b"acgt"), b"acgu");
    }

    #[test]
    fn rna_no_t() {
        assert_eq!(dna_to_rna(b"ACGA"), b"ACGA");
    }

    #[test]
    fn rna_passthrough_non_dna() {
        assert_eq!(dna_to_rna(b"AFCGA"), b"AFCGA");
    }

    //Reverse Complement
    #[test]
    fn reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn reverse_complement_single_bases() {
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"T"), b"A");
        assert_eq!(reverse_complement(b"G"), b"C");
        assert_eq!(reverse_complement(b"C"), b"G");
        assert_eq!(reverse_complement(b"U"), b"A");
    }

    #[test]
    fn reverse_complement_ambiguous_passthrough() {
        assert_eq!(reverse_complement(b"N"), b"N");
        assert_eq!(reverse_complement(b"n"), b"n");
    }

    #[test]
    fn reverse_complement_idempotent() {
        let seq = b"ACGTACGT";
        assert_eq!(reverse_complement(&reverse_complement(seq)), seq);
    }

    #[test]
    fn reverse_complement_pure_dna() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
    }

    #[test]
    fn reverse() {
        assert_eq!(reverse_complement(b"AGCTUagctuNn"), b"nNaagctAAGCT")
    }

    // LengthFilter tests
    #[test]
    fn length_filter_gt() {
        let f = LengthFilter {
            op: LengthOp::Gt,
            value: 10,
        };
        assert!(f.matches(11));
        assert!(!f.matches(10));
        assert!(!f.matches(9));
    }

    #[test]
    fn length_filter_eq() {
        let f = LengthFilter {
            op: LengthOp::Eq,
            value: 10,
        };
        assert!(f.matches(10));
        assert!(!f.matches(11));
        assert!(!f.matches(9));
    }

    // SequenceFilter tests
    #[test]
    fn sequence_filter_contains() {
        let f = SequenceFilter {
            op: SequenceOp::Contains,
            pattern: "ACGT".to_string(),
        };
        assert!(f.like(b"GGACGTGG"));
        assert!(!f.like(b"GGAAGGGG"));
    }

    #[test]
    fn sequence_filter_starts_with() {
        let f = SequenceFilter {
            op: SequenceOp::StartsWith,
            pattern: "ACGT".to_string(),
        };
        assert!(f.like(b"ACGTGGGG"));
        assert!(!f.like(b"GGACGTGG"));
    }

    #[test]
    fn sequence_filter_ends_with() {
        let f = SequenceFilter {
            op: SequenceOp::EndsWith,
            pattern: "ACGT".to_string(),
        };
        assert!(f.like(b"GGGGACGT"));
        assert!(!f.like(b"GGACGTGG"));
    }

    // is_valid_dna
    #[test]
    fn valid_dna_basic() {
        assert!(is_valid_dna(b"ACGT"));
    }

    #[test]
    fn valid_dna_lowercase() {
        assert!(is_valid_dna(b"acgt"));
    }

    #[test]
    fn valid_dna_mixed_case() {
        assert!(is_valid_dna(b"AcGt"));
    }

    #[test]
    fn valid_dna_with_n() {
        assert!(is_valid_dna(b"ACGTN"));
        assert!(is_valid_dna(b"acgtn"));
    }

    #[test]
    fn valid_dna_empty() {
        assert!(is_valid_dna(b""));
    }

    #[test]
    fn valid_dna_rejects_u() {
        assert!(!is_valid_dna(b"ACGTU"));
    }

    #[test]
    fn valid_dna_rejects_invalid() {
        assert!(!is_valid_dna(b"ACGTZ"));
        assert!(!is_valid_dna(b"ACGT1"));
        assert!(!is_valid_dna(b"ACGT!"));
    }

    #[test]
    fn valid_dna_rejects_on_first_invalid() {
        // invalid base at start — should short circuit
        assert!(!is_valid_dna(b"ZACGT"));
    }

    #[test]
    fn valid_dna_single_base() {
        assert!(is_valid_dna(b"A"));
        assert!(is_valid_dna(b"C"));
        assert!(is_valid_dna(b"G"));
        assert!(is_valid_dna(b"T"));
        assert!(is_valid_dna(b"N"));
        assert!(!is_valid_dna(b"U"));
        assert!(!is_valid_dna(b"Z"));
    }

    // is_valid_rna
    #[test]
    fn valid_rna_basic() {
        assert!(is_valid_rna(b"ACGU"));
    }

    #[test]
    fn valid_rna_lowercase() {
        assert!(is_valid_rna(b"acgu"));
    }

    #[test]
    fn valid_rna_mixed_case() {
        assert!(is_valid_rna(b"AcGu"));
    }

    #[test]
    fn valid_rna_with_n() {
        assert!(is_valid_rna(b"ACGUN"));
        assert!(is_valid_rna(b"acgun"));
    }

    #[test]
    fn valid_rna_empty() {
        assert!(is_valid_rna(b""));
    }

    #[test]
    fn valid_rna_rejects_t() {
        assert!(!is_valid_rna(b"ACGT"));
    }

    #[test]
    fn valid_rna_rejects_invalid() {
        assert!(!is_valid_rna(b"ACGUZ"));
        assert!(!is_valid_rna(b"ACGU1"));
        assert!(!is_valid_rna(b"ACGU!"));
    }

    #[test]
    fn valid_rna_single_base() {
        assert!(is_valid_rna(b"A"));
        assert!(is_valid_rna(b"C"));
        assert!(is_valid_rna(b"G"));
        assert!(is_valid_rna(b"U"));
        assert!(is_valid_rna(b"N"));
        assert!(!is_valid_rna(b"T"));
        assert!(!is_valid_rna(b"Z"));
    }
}
