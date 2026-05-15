use std::fs::File;
use std::io::{BufWriter, Write};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let count: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(50_000);
    let min_len: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(300);
    let max_len: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(600);

    let fa_path = "bench_reads.fa";
    let fai_path = "bench_reads.fa.fai";
    let nofai_path = "bench_reads_nofai.fa";

    let mut fa = BufWriter::new(File::create(fa_path).unwrap());
    let mut fai = BufWriter::new(File::create(fai_path).unwrap());

    let mut seed: u64 = 99991;
    let mut rand = move || -> u64 {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };

    let bases = b"ACGT";
    let mut offset: u64 = 0;

    for i in 0..count {
        let length = (rand() % (max_len - min_len + 1) as u64) as usize + min_len;
        let header = format!(">seq{i}\n");
        let seq_offset = offset + header.len() as u64;

        fa.write_all(header.as_bytes()).unwrap();
        let seq: Vec<u8> = (0..length).map(|_| bases[(rand() % 4) as usize]).collect();
        fa.write_all(&seq).unwrap();
        fa.write_all(b"\n").unwrap();

        writeln!(fai, "seq{i}\t{length}\t{seq_offset}\t{length}\t{}", length + 1).unwrap();

        offset = seq_offset + length as u64 + 1;
    }

    fa.flush().unwrap();
    fai.flush().unwrap();

    std::fs::copy(fa_path, nofai_path).unwrap();
    eprintln!("Generated {count} records ({min_len}–{max_len} bp) → {fa_path} + {fai_path}");
    eprintln!("Copied to {nofai_path} (no FAI index)");
}
