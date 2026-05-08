use sqlite3_ext::*;

use crate::filters::ExecPlan;

pub trait SequenceRecord: Clone {
    fn identifier_bytes(&self) -> &[u8];
    fn description_bytes(&self) -> Option<&[u8]>;
    fn sequence_bytes(&self) -> &[u8];
}
pub trait SequenceReader {
    type Record: SequenceRecord;
    fn next(&mut self) -> Option<Result<Self::Record>>;
}

pub struct SequenceCursor<R: SequenceReader> {
    pub plan: ExecPlan,
    pub fallback_filename: Option<String>,
    pub reader: Option<R>,
    pub current: Option<R::Record>,
    pub rowid: i64,
    pub done: bool,
}
