use sqlite3_ext::{function::FunctionOptions, *};
pub mod functions;
use functions::*;
pub mod filters;
mod reader;
use reader::*;
mod fasta;
use fasta::*;
mod fastq;
use fastq::*;
pub mod fai;

#[sqlite3_ext_main]
pub fn init(db: &Connection) -> Result<()> {
    db.create_module("fasta", FastaModule::module(), ())?;
    db.create_module("fastq", FastqModule::module(), ())?;
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
                .ok_or("base_count requires a non-empty base argument")?;
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
        "complement",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let seq = complement(seq.as_bytes());
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
    db.create_scalar_function(
        "min_quality",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let quality = args[0].get_str()?;
            let valid = min_quality(quality.as_bytes());
            ctx.set_result(valid)
        },
    )?;
    db.create_scalar_function(
        "max_quality",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let quality = args[0].get_str()?;
            let valid = max_quality(quality.as_bytes());
            ctx.set_result(valid)
        },
    )?;
    db.create_scalar_function(
        "mean_quality",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let quality = args[0].get_str()?;
            let valid = mean_quality(quality.as_bytes());
            ctx.set_result(valid)
        },
    )?;
    db.create_scalar_function(
        "base_composition",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let comp = base_composition(seq.as_bytes());
            ctx.set_result(comp.to_string())
        },
    )?;
    db.create_scalar_function(
        "has_stop_codon",
        &FunctionOptions::default().set_n_args(1),
        |ctx, args| {
            let seq = args[0].get_str()?;
            let result = has_stop_codon(seq.as_bytes());
            ctx.set_result(result)
        },
    )?;
    db.create_aggregate_function::<(), N50Accumulator>(
        "n50",
        &FunctionOptions::default().set_n_args(1),
        (),
    )?;
    db.create_aggregate_function::<(), BaseTotalsAggregator>(
        "base_totals",
        &FunctionOptions::default().set_n_args(1),
        (),
    )?;
    Ok(())
}
