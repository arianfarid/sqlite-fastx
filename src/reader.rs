use sqlite3_ext::*;

use crate::filters::ExecPlan;

pub trait SequenceRecord: Clone {
    fn identifier_bytes(&self) -> &[u8];
    fn description_bytes(&self) -> Option<&[u8]>;
    fn sequence_bytes(&self) -> &[u8];
    fn quality_bytes(&self) -> Option<&[u8]>;
}
pub trait SequenceReader {
    type Record: SequenceRecord;
    fn next(&mut self) -> Option<Result<Self::Record>>;
    fn lookup_offset(_fai_path: &str, _id: &str) -> Option<u64> {
        None
    }
}

pub struct SequenceCursor<R: SequenceReader> {
    pub plan: ExecPlan,
    pub fallback_filename: Option<String>,
    pub fai_path: Option<String>,
    pub reader: Option<R>,
    pub current: Option<R::Record>,
    pub rowid: i64,
    pub done: bool,
    pub exit_early: bool,
}

impl<R: SequenceReader> SequenceCursor<R> {
    pub fn determine_strategy(
        &mut self,
        index_str: Option<&str>,
        args: &mut [&mut ValueRef],
    ) -> Result<ReadStrategy> {
        let Some(descriptor) = index_str else {
            return Ok(ReadStrategy::Stream);
        };

        let id_arg_idx = descriptor.split(',').position(|t| t == "id:Eq");
        let Some(idx) = id_arg_idx else {
            return Ok(ReadStrategy::Stream);
        };

        let Some(ref fai_path) = self.fai_path else {
            return Ok(ReadStrategy::Stream);
        };
        let id = args[idx].get_str()?.to_string();

        match R::lookup_offset(fai_path, &id) {
            Some(offset) => Ok(ReadStrategy::SeekToOffset(offset)),
            None => Ok(ReadStrategy::Stream),
        }
    }
}
pub enum ReadStrategy {
    Stream,
    SeekToOffset(u64),
}
