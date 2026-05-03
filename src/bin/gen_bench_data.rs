use std::io::{self, Write};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let count: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(10_000);
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    // Simple LCG random number generator — no dependencies needed
    let mut seed: u64 = 12345;
    let mut rand = move || -> u64 {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };

    let bases = b"ACGT";

    for i in 0..count {
        let length = (rand() % 451 + 50) as usize; // 50–500
        write!(out, ">seq{}\n", i).unwrap();
        for _ in 0..length {
            out.write_all(&[bases[(rand() % 4) as usize]]).unwrap();
        }
        writeln!(out).unwrap();
    }
}
