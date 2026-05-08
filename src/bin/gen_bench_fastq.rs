// src/bin/gen_bench_fastq.rs
use std::io::{self, Write};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let count: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(10_000);
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    // Simple LCG random number generator
    let mut seed: u64 = 98765;
    let mut rand = move || -> u64 {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };

    let bases = b"ACGT";
    // Phred+33 quality scores — printable ASCII from '!' (33) to 'J' (74)
    // representing quality scores 0-41
    let quals = b"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ";

    for i in 0..count {
        let length = (rand() % 151 + 50) as usize; // 50-200 typical short read length

        let seq: Vec<u8> = (0..length).map(|_| bases[(rand() % 4) as usize]).collect();

        let qual: Vec<u8> = (0..length)
            .map(|_| quals[(rand() % quals.len() as u64) as usize])
            .collect();

        writeln!(out, "@seq{}", i).unwrap();
        out.write_all(&seq).unwrap();
        writeln!(out).unwrap();
        writeln!(out, "+").unwrap();
        out.write_all(&qual).unwrap();
        writeln!(out).unwrap();
    }
}
