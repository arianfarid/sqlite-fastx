use sqlite3_ext::*;

use crate::{functions::compute_gc, reader::SequenceRecord};

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
    pub fn matches(&self, len: i64) -> bool {
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
    pub fn matches(&self, gc: f64) -> bool {
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
impl SequenceOp {
    // pub fn as_str(&self) -> &'static str {
    //     match self {
    //         SequenceOp::Contains => "Contains",
    //         SequenceOp::StartsWith => "StartsWith",
    //         SequenceOp::EndsWith => "EndsWith",
    //         SequenceOp::Eq => "Eq",
    //     }
    // }
}
pub struct SequenceFilter {
    pub op: SequenceOp,
    pub pattern: String,
}
impl SequenceFilter {
    pub fn like(&self, seq: &[u8]) -> bool {
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

enum Predicate {
    Length(LengthFilter),
    SequenceLike(SequenceFilter),
    GC(GCFilter),
    // TODO: Substring
}
impl Predicate {
    fn matches<S: SequenceRecord>(&self, record: &S) -> bool {
        match self {
            Predicate::Length(f) => f.matches(record.sequence_bytes().len() as i64),
            Predicate::GC(f) => f.matches(compute_gc(record.sequence_bytes())),
            Predicate::SequenceLike(s) => s.like(record.sequence_bytes()),
        }
    }
}

pub struct ExecPlan {
    predicates: Vec<Predicate>,
}
impl ExecPlan {
    pub fn new() -> ExecPlan {
        ExecPlan { predicates: vec![] }
    }
    pub fn matches<S: SequenceRecord>(&self, record: &S) -> bool {
        for pred in &self.predicates {
            if !pred.matches(record) {
                return false;
            }
        }
        true
    }
}

pub fn parse_plan(index_str: Option<&str>, args: &mut [&mut ValueRef]) -> Result<ExecPlan> {
    let Some(descriptor) = index_str else {
        return Ok(ExecPlan { predicates: vec![] });
    };

    let mut predicates = vec![];
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
            ("sequence", "Like") => {
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
        assert!(f.matches(11));
        assert!(!f.matches(10));
        assert!(!f.matches(9));
    }

    #[test]
    fn length_filter_eq() {
        let f = LengthFilter {
            op: CompareOp::Eq,
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

    // GC Filter
    #[test]
    fn gc_filter_gt() {
        let f = GCFilter {
            op: CompareOp::Gt,
            value: 0.5,
        };
        assert!(f.matches(0.6));
        assert!(!f.matches(0.5));
        assert!(!f.matches(0.4));
    }

    #[test]
    fn gc_filter_ge() {
        let f = GCFilter {
            op: CompareOp::Ge,
            value: 0.5,
        };
        assert!(f.matches(0.6));
        assert!(f.matches(0.5));
        assert!(!f.matches(0.4));
    }

    #[test]
    fn gc_filter_lt() {
        let f = GCFilter {
            op: CompareOp::Lt,
            value: 0.5,
        };
        assert!(f.matches(0.4));
        assert!(!f.matches(0.5));
        assert!(!f.matches(0.6));
    }

    #[test]
    fn gc_filter_le() {
        let f = GCFilter {
            op: CompareOp::Le,
            value: 0.5,
        };
        assert!(f.matches(0.4));
        assert!(f.matches(0.5));
        assert!(!f.matches(0.6));
    }

    #[test]
    fn gc_filter_eq() {
        let f = GCFilter {
            op: CompareOp::Eq,
            value: 0.5,
        };
        assert!(f.matches(0.5));
        assert!(!f.matches(0.5 + f64::EPSILON * 2.0));
        assert!(!f.matches(0.4));
    }

    #[test]
    fn gc_filter_boundaries() {
        let f = GCFilter {
            op: CompareOp::Ge,
            value: 0.0,
        };
        assert!(f.matches(0.0));
        assert!(f.matches(1.0));

        let f = GCFilter {
            op: CompareOp::Le,
            value: 1.0,
        };
        assert!(f.matches(0.0));
        assert!(f.matches(1.0));
    }

    #[test]
    fn gc_filter_pure_gc() {
        // 1.0 GC content
        let f = GCFilter {
            op: CompareOp::Eq,
            value: 1.0,
        };
        assert!(f.matches(1.0));
        assert!(!f.matches(0.99));
    }

    #[test]
    fn gc_filter_no_gc() {
        // 0.0 GC content
        let f = GCFilter {
            op: CompareOp::Eq,
            value: 0.0,
        };
        assert!(f.matches(0.0));
        assert!(!f.matches(0.01));
    }
}
