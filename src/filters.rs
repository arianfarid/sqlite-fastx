use sqlite3_ext::*;

use crate::reader::SequenceRecord;

#[repr(i32)]
pub enum LengthOp {
    Gt,
    Ge,
    Lt,
    Le,
    Eq,
}

impl LengthOp {
    pub fn as_str(&self) -> &'static str {
        match self {
            LengthOp::Gt => "Gt",
            LengthOp::Ge => "Ge",
            LengthOp::Lt => "Lt",
            LengthOp::Le => "Le",
            LengthOp::Eq => "Eq",
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

pub struct LengthFilter {
    pub op: LengthOp,
    pub value: i64,
}
impl LengthFilter {
    pub fn matches(&self, len: i64) -> bool {
        match self.op {
            LengthOp::Gt => len > self.value,
            LengthOp::Ge => len >= self.value,
            LengthOp::Lt => len < self.value,
            LengthOp::Le => len <= self.value,
            LengthOp::Eq => len == self.value,
        }
    }
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

#[cfg(test)]
mod tests {
    use super::*;
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
}
