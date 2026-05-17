use sqlite3_ext::*;

use crate::{
    functions::{compute_gc, mean_quality, min_quality},
    reader::SequenceRecord,
};

#[repr(i32)]
pub enum CompareOp {
    Gt,
    Ge,
    Lt,
    Le,
    Eq,
}

impl CompareOp {
    pub fn as_str(&self) -> &'static str {
        match self {
            CompareOp::Gt => "Gt",
            CompareOp::Ge => "Ge",
            CompareOp::Lt => "Lt",
            CompareOp::Le => "Le",
            CompareOp::Eq => "Eq",
        }
    }
}
pub struct LengthFilter {
    pub op: CompareOp,
    pub value: i64,
}
impl LengthFilter {
    pub fn eval(&self, len: i64) -> bool {
        match self.op {
            CompareOp::Gt => len > self.value,
            CompareOp::Ge => len >= self.value,
            CompareOp::Lt => len < self.value,
            CompareOp::Le => len <= self.value,
            CompareOp::Eq => len == self.value,
        }
    }
}
pub struct GCFilter {
    pub op: CompareOp,
    pub value: f64,
}
impl GCFilter {
    pub fn eval(&self, gc: f64) -> bool {
        match self.op {
            CompareOp::Gt => gc > self.value,
            CompareOp::Ge => gc >= self.value,
            CompareOp::Lt => gc < self.value,
            CompareOp::Le => gc <= self.value,
            CompareOp::Eq => (gc - self.value).abs() < f64::EPSILON,
        }
    }
}
pub enum SequenceOp {
    Contains,
    StartsWith,
    EndsWith,
    Eq,
}

pub struct EqFilter {
    pub op: SequenceOp,
    pub pattern: String,
}
impl EqFilter {
    pub fn eval(&self, val: &[u8]) -> bool {
        val == self.pattern.as_bytes()
    }
}
pub struct LikeFilter {
    pub op: SequenceOp,
    pub pattern: String,
}
impl LikeFilter {
    pub fn eval(&self, seq: &[u8]) -> bool {
        let val = seq.to_ascii_uppercase();
        match self.op {
            SequenceOp::Contains => memchr::memmem::find(&val, self.pattern.as_bytes()).is_some(),
            SequenceOp::StartsWith => val.starts_with(self.pattern.as_bytes()),
            SequenceOp::EndsWith => val.ends_with(self.pattern.as_bytes()),
            SequenceOp::Eq => val == self.pattern.as_bytes(),
        }
    }
}

enum Predicate {
    IDLike(LikeFilter),
    IDEq(EqFilter),
    Description(LikeFilter),
    Length(LengthFilter),
    Sequence(LikeFilter),
    GC(GCFilter),
    MinQuality(LengthFilter),
    MeanQuality(GCFilter),
}
impl Predicate {
    fn eval<S: SequenceRecord>(&self, record: &S) -> bool {
        match self {
            Predicate::IDEq(f) => f.eval(record.identifier_bytes()),
            Predicate::IDLike(f) => f.eval(record.identifier_bytes()),
            Predicate::Description(f) => record
                .description_bytes()
                .map(|d| f.eval(d))
                .unwrap_or(false),
            Predicate::Length(f) => f.eval(record.sequence_bytes().len() as i64),
            Predicate::GC(f) => f.eval(compute_gc(record.sequence_bytes())),
            Predicate::Sequence(s) => s.eval(record.sequence_bytes()),
            Predicate::MeanQuality(f) => {
                f.eval(mean_quality(record.quality_bytes().unwrap_or_default()))
            }
            Predicate::MinQuality(f) => {
                f.eval(min_quality(record.quality_bytes().unwrap_or_default()))
            }
        }
    }
}

pub struct ExecPlan {
    predicates: Vec<Predicate>,
    pub unique: bool,
}
impl Default for ExecPlan {
    fn default() -> Self {
        Self::new()
    }
}
impl ExecPlan {
    pub fn new() -> ExecPlan {
        ExecPlan {
            predicates: vec![],
            unique: false,
        }
    }
    pub fn eval<S: SequenceRecord>(&self, record: &S) -> bool {
        for pred in &self.predicates {
            if !pred.eval(record) {
                return false;
            }
        }
        true
    }
}

pub fn parse_plan(index_str: Option<&str>, args: &mut [&mut ValueRef]) -> Result<ExecPlan> {
    let Some(descriptor) = index_str else {
        return Ok(ExecPlan {
            predicates: vec![],
            unique: false,
        });
    };

    let mut predicates = vec![];
    let mut unique = false;
    for (i, token) in descriptor.split(',').enumerate() {
        if i >= args.len() {
            break;
        }
        let arg = &mut args[i];

        let mut parts = token.splitn(2, ":");
        let col = parts.next().unwrap_or("");
        let op = parts.next().unwrap_or("");

        match (col, op) {
            ("length", "Gt") => predicates.push(Predicate::Length(LengthFilter {
                op: CompareOp::Gt,
                value: arg.get_i64(),
            })),
            ("length", "Ge") => predicates.push(Predicate::Length(LengthFilter {
                op: CompareOp::Ge,
                value: arg.get_i64(),
            })),
            ("length", "Lt") => predicates.push(Predicate::Length(LengthFilter {
                op: CompareOp::Lt,
                value: arg.get_i64(),
            })),
            ("length", "Le") => predicates.push(Predicate::Length(LengthFilter {
                op: CompareOp::Le,
                value: arg.get_i64(),
            })),
            ("length", "Eq") => predicates.push(Predicate::Length(LengthFilter {
                op: CompareOp::Eq,
                value: arg.get_i64(),
            })),
            ("gc_content", "Gt") => predicates.push(Predicate::GC(GCFilter {
                op: CompareOp::Gt,
                value: arg.get_f64(),
            })),
            ("gc_content", "Ge") => predicates.push(Predicate::GC(GCFilter {
                op: CompareOp::Ge,
                value: arg.get_f64(),
            })),
            ("gc_content", "Lt") => predicates.push(Predicate::GC(GCFilter {
                op: CompareOp::Lt,
                value: arg.get_f64(),
            })),
            ("gc_content", "Le") => predicates.push(Predicate::GC(GCFilter {
                op: CompareOp::Le,
                value: arg.get_f64(),
            })),
            ("gc_content", "Eq") => predicates.push(Predicate::GC(GCFilter {
                op: CompareOp::Eq,
                value: arg.get_f64(),
            })),
            ("id", "Eq") => {
                let raw = arg.get_str()?.to_string();
                // Per fasta rules, id = 'some_id' must be unique
                unique = true;
                predicates.push(Predicate::IDEq(EqFilter {
                    op: SequenceOp::Eq,
                    pattern: raw,
                }))
            }
            ("id", "Like") => {
                let raw = arg.get_str()?.to_string();
                let (op, pattern) = parse_like_pattern(&raw);
                predicates.push(Predicate::IDLike(LikeFilter { op, pattern }))
            }
            ("description", "Like") => {
                let raw = arg.get_str()?.to_string();
                let (op, pattern) = parse_like_pattern(&raw);
                predicates.push(Predicate::Description(LikeFilter { op, pattern }))
            }
            ("sequence", "Like") => {
                let raw = arg.get_str()?.to_string();
                let (op, pattern) = parse_like_pattern(&raw);
                predicates.push(Predicate::Sequence(LikeFilter { op, pattern }))
            }
            ("mean_quality", "Gt") => predicates.push(Predicate::MeanQuality(GCFilter {
                op: CompareOp::Gt,
                value: arg.get_f64(),
            })),
            ("mean_quality", "Ge") => predicates.push(Predicate::MeanQuality(GCFilter {
                op: CompareOp::Ge,
                value: arg.get_f64(),
            })),
            ("mean_quality", "Lt") => predicates.push(Predicate::MeanQuality(GCFilter {
                op: CompareOp::Lt,
                value: arg.get_f64(),
            })),
            ("mean_quality", "Le") => predicates.push(Predicate::MeanQuality(GCFilter {
                op: CompareOp::Le,
                value: arg.get_f64(),
            })),
            ("mean_quality", "Eq") => predicates.push(Predicate::MeanQuality(GCFilter {
                op: CompareOp::Eq,
                value: arg.get_f64(),
            })),
            ("min_quality", "Gt") => predicates.push(Predicate::MinQuality(LengthFilter {
                op: CompareOp::Gt,
                value: arg.get_i64(),
            })),
            ("min_quality", "Ge") => predicates.push(Predicate::MinQuality(LengthFilter {
                op: CompareOp::Ge,
                value: arg.get_i64(),
            })),
            ("min_quality", "Lt") => predicates.push(Predicate::MinQuality(LengthFilter {
                op: CompareOp::Lt,
                value: arg.get_i64(),
            })),
            ("min_quality", "Le") => predicates.push(Predicate::MinQuality(LengthFilter {
                op: CompareOp::Le,
                value: arg.get_i64(),
            })),
            ("min_quality", "Eq") => predicates.push(Predicate::MinQuality(LengthFilter {
                op: CompareOp::Eq,
                value: arg.get_i64(),
            })),
            _ => {}
        };
    }
    Ok(ExecPlan { predicates, unique })
}

fn parse_like_pattern(raw: &str) -> (SequenceOp, String) {
    let starts_with_wild = raw.starts_with('%');
    let ends_with_wild = raw.ends_with('%');
    let pattern = raw.trim_matches('%').to_ascii_uppercase();
    let op = match (starts_with_wild, ends_with_wild) {
        (true, true) => SequenceOp::Contains,
        (true, false) => SequenceOp::EndsWith,
        (false, true) => SequenceOp::StartsWith,
        (false, false) => SequenceOp::Eq,
    };
    (op, pattern)
}

#[cfg(test)]
mod tests {
    use super::*;
    // LengthFilter tests
    #[test]
    fn length_filter_gt() {
        let f = LengthFilter {
            op: CompareOp::Gt,
            value: 10,
        };
        assert!(f.eval(11));
        assert!(!f.eval(10));
        assert!(!f.eval(9));
    }

    #[test]
    fn length_filter_eq() {
        let f = LengthFilter {
            op: CompareOp::Eq,
            value: 10,
        };
        assert!(f.eval(10));
        assert!(!f.eval(11));
        assert!(!f.eval(9));
    }

    // LikeFilter tests
    #[test]
    fn like_filter_contains() {
        let f = LikeFilter {
            op: SequenceOp::Contains,
            pattern: "ACGT".to_string(),
        };
        assert!(f.eval(b"GGACGTGG"));
        assert!(!f.eval(b"GGAAGGGG"));
    }

    #[test]
    fn like_filter_starts_with() {
        let f = LikeFilter {
            op: SequenceOp::StartsWith,
            pattern: "ACGT".to_string(),
        };
        assert!(f.eval(b"ACGTGGGG"));
        assert!(!f.eval(b"GGACGTGG"));
    }

    #[test]
    fn like_filter_ends_with() {
        let f = LikeFilter {
            op: SequenceOp::EndsWith,
            pattern: "ACGT".to_string(),
        };
        assert!(f.eval(b"GGGGACGT"));
        assert!(!f.eval(b"GGACGTGG"));
    }

    #[test]
    fn like_filter_case_insensitive_contains() {
        let f = LikeFilter {
            op: SequenceOp::Contains,
            pattern: "ACGT".to_string(),
        };
        assert!(f.eval(b"ggacgtgg"));
        assert!(f.eval(b"ggACGTgg"));
    }

    #[test]
    fn like_filter_case_insensitive_starts_with() {
        let f = LikeFilter {
            op: SequenceOp::StartsWith,
            pattern: "ACGT".to_string(),
        };
        assert!(f.eval(b"acgtgggg"));
    }

    #[test]
    fn like_filter_case_insensitive_ends_with() {
        let f = LikeFilter {
            op: SequenceOp::EndsWith,
            pattern: "ACGT".to_string(),
        };
        assert!(f.eval(b"ggggacgt"));
    }

    #[test]
    fn like_filter_case_insensitive_eq() {
        let f = LikeFilter {
            op: SequenceOp::Eq,
            pattern: "ACGT".to_string(),
        };
        assert!(f.eval(b"acgt"));
        assert!(f.eval(b"ACGT"));
        assert!(f.eval(b"AcGt"));
    }

    // GC Filter
    #[test]
    fn gc_filter_gt() {
        let f = GCFilter {
            op: CompareOp::Gt,
            value: 0.5,
        };
        assert!(f.eval(0.6));
        assert!(!f.eval(0.5));
        assert!(!f.eval(0.4));
    }

    #[test]
    fn gc_filter_ge() {
        let f = GCFilter {
            op: CompareOp::Ge,
            value: 0.5,
        };
        assert!(f.eval(0.6));
        assert!(f.eval(0.5));
        assert!(!f.eval(0.4));
    }

    #[test]
    fn gc_filter_lt() {
        let f = GCFilter {
            op: CompareOp::Lt,
            value: 0.5,
        };
        assert!(f.eval(0.4));
        assert!(!f.eval(0.5));
        assert!(!f.eval(0.6));
    }

    #[test]
    fn gc_filter_le() {
        let f = GCFilter {
            op: CompareOp::Le,
            value: 0.5,
        };
        assert!(f.eval(0.4));
        assert!(f.eval(0.5));
        assert!(!f.eval(0.6));
    }

    #[test]
    fn gc_filter_eq() {
        let f = GCFilter {
            op: CompareOp::Eq,
            value: 0.5,
        };
        assert!(f.eval(0.5));
        assert!(!f.eval(0.5 + f64::EPSILON * 2.0));
        assert!(!f.eval(0.4));
    }

    #[test]
    fn gc_filter_boundaries() {
        let f = GCFilter {
            op: CompareOp::Ge,
            value: 0.0,
        };
        assert!(f.eval(0.0));
        assert!(f.eval(1.0));

        let f = GCFilter {
            op: CompareOp::Le,
            value: 1.0,
        };
        assert!(f.eval(0.0));
        assert!(f.eval(1.0));
    }

    #[test]
    fn gc_filter_pure_gc() {
        // 1.0 GC content
        let f = GCFilter {
            op: CompareOp::Eq,
            value: 1.0,
        };
        assert!(f.eval(1.0));
        assert!(!f.eval(0.99));
    }

    #[test]
    fn gc_filter_no_gc() {
        // 0.0 GC content
        let f = GCFilter {
            op: CompareOp::Eq,
            value: 0.0,
        };
        assert!(f.eval(0.0));
        assert!(!f.eval(0.01));
    }

    // parse_like_pattern
    #[test]
    fn parse_like_both_wildcards() {
        let (op, pattern) = parse_like_pattern("%ACGT%");
        assert!(matches!(op, SequenceOp::Contains));
        assert_eq!(pattern, "ACGT");
    }

    #[test]
    fn parse_like_leading_wildcard() {
        let (op, pattern) = parse_like_pattern("%ACGT");
        assert!(matches!(op, SequenceOp::EndsWith));
        assert_eq!(pattern, "ACGT");
    }

    #[test]
    fn parse_like_trailing_wildcard() {
        let (op, pattern) = parse_like_pattern("ACGT%");
        assert!(matches!(op, SequenceOp::StartsWith));
        assert_eq!(pattern, "ACGT");
    }

    #[test]
    fn parse_like_no_wildcard() {
        let (op, pattern) = parse_like_pattern("ACGT");
        assert!(matches!(op, SequenceOp::Eq));
        assert_eq!(pattern, "ACGT");
    }

    #[test]
    fn parse_like_uppercases_pattern() {
        let (_, pattern) = parse_like_pattern("%acgt%");
        assert_eq!(pattern, "ACGT");
    }
}
