use criterion::{Criterion, criterion_group, criterion_main};
use sqlite_fastx::filters::{LikeFilter, SequenceOp};
use sqlite_fastx::init;
use sqlite3_ext::{Database, FallibleIteratorMut, FromValue};
use std::hint::black_box;

const BENCH_FA: &str = "bench.fa";
const BENCH_WIDE_FAI: &str = "bench_wide.fa";
const BENCH_WIDE_NOFAI: &str = "bench_wide_nofai.fa";
const WIDE_COUNT: usize = 20;

fn setup() -> Database {
    let db = Database::open(":memory:").unwrap();
    init(&db).unwrap();
    db.execute(
        &format!("CREATE VIRTUAL TABLE fa USING fasta('{}')", BENCH_FA),
        (),
    )
    .unwrap();
    db
}

fn count(db: &Database, sql: &str) -> i64 {
    let mut stmt = db.prepare(sql).unwrap();
    let rows = stmt.query(()).unwrap();
    let row = rows.next().unwrap().unwrap();
    row[0].get_i64()
}

// --- SQL-level pushdown benchmarks ---

fn bench_full_scan(c: &mut Criterion) {
    let db = setup();
    c.bench_function("full_scan_10k", |b| {
        b.iter(|| black_box(count(&db, "SELECT COUNT(*) FROM fa")))
    });
}

fn bench_length_pushdown(c: &mut Criterion) {
    let db = setup();
    let mut group = c.benchmark_group("length_gt_100");
    group.bench_function("pushdown", |b| {
        b.iter(|| black_box(count(&db, "SELECT COUNT(*) FROM fa WHERE length > 100")))
    });
    group.bench_function("no_pushdown", |b| {
        b.iter(|| black_box(count(&db, "SELECT COUNT(*) FROM fa WHERE length + 0 > 100")))
    });
    group.finish();
}

fn bench_length_range(c: &mut Criterion) {
    let db = setup();
    let mut group = c.benchmark_group("length_gt_100_lt_200");
    group.bench_function("pushdown", |b| {
        b.iter(|| {
            black_box(count(
                &db,
                "SELECT COUNT(*) FROM fa WHERE length > 100 AND length < 200",
            ))
        })
    });
    group.bench_function("no_pushdown", |b| {
        b.iter(|| {
            black_box(count(
                &db,
                "SELECT COUNT(*) FROM fa WHERE length + 0 > 100 AND length + 0 < 200",
            ))
        })
    });
    group.finish();
}

fn bench_sequence_contains(c: &mut Criterion) {
    let db = setup();
    let mut group = c.benchmark_group("sequence_contains_ACGT");
    group.bench_function("pushdown", |b| {
        b.iter(|| {
            black_box(count(
                &db,
                "SELECT COUNT(*) FROM fa WHERE sequence LIKE '%ACGT%'",
            ))
        })
    });
    group.bench_function("no_pushdown", |b| {
        b.iter(|| {
            black_box(count(
                &db,
                "SELECT COUNT(*) FROM fa WHERE CAST(sequence AS TEXT) LIKE '%ACGT%'",
            ))
        })
    });
    group.finish();
}

fn bench_combined(c: &mut Criterion) {
    let db = setup();
    let mut group = c.benchmark_group("length_and_sequence");
    group.bench_function("pushdown", |b| {
        b.iter(|| {
            black_box(count(
                &db,
                "SELECT COUNT(*) FROM fa WHERE length > 100 AND sequence LIKE '%ACGT%'",
            ))
        })
    });
    group.bench_function("no_pushdown", |b| {
        b.iter(|| {
            black_box(count(
                &db,
                "SELECT COUNT(*) FROM fa WHERE length + 0 > 100 AND CAST(sequence AS TEXT) LIKE '%ACGT%'",
            ))
        })
    });
    group.finish();
}

// --- Raw function benchmarks ---

fn bench_sequence_search(c: &mut Criterion) {
    let short_no_match: Vec<u8> = b"A".repeat(150);
    let short_early_match: Vec<u8> = {
        let mut v = b"ACGT".to_vec();
        v.extend(b"A".repeat(146));
        v
    };
    let short_late_match: Vec<u8> = {
        let mut v = b"A".repeat(146);
        v.extend_from_slice(b"ACGT");
        v
    };

    let long_no_match: Vec<u8> = b"A".repeat(10_000);
    let long_early_match: Vec<u8> = {
        let mut v = b"ACGT".to_vec();
        v.extend(b"A".repeat(9_996));
        v
    };
    let long_late_match: Vec<u8> = {
        let mut v = b"A".repeat(9_996);
        v.extend_from_slice(b"ACGT");
        v
    };
    let long_multiple_matches: Vec<u8> = {
        let mut v = vec![];
        for _ in 0..100 {
            v.extend(b"A".repeat(96));
            v.extend_from_slice(b"ACGT");
        }
        v
    };
    let long_realistic: Vec<u8> = {
        let bases = b"GCTAGCTAGCTAGC";
        let mut v: Vec<u8> = (0..5_000).map(|i| bases[i % bases.len()]).collect();
        v.extend_from_slice(b"ACGT");
        v.extend((0..4_996).map(|i| bases[i % bases.len()]));
        v
    };
    let long_near_miss: Vec<u8> = {
        let mut v: Vec<u8> = b"ACG".repeat(3_333);
        v.push(b'A');
        v
    };

    let vlong_no_match: Vec<u8> = b"A".repeat(100_000);
    let vlong_late_match: Vec<u8> = {
        let mut v = b"A".repeat(99_996);
        v.extend_from_slice(b"ACGT");
        v
    };
    let vlong_near_miss: Vec<u8> = {
        let mut v: Vec<u8> = b"ACG".repeat(33_333);
        v.push(b'A');
        v
    };
    let vlong_realistic: Vec<u8> = {
        let bases = b"GCTAGCTAGCTAGC";
        let mut v: Vec<u8> = (0..50_000).map(|i| bases[i % bases.len()]).collect();
        v.extend_from_slice(b"ACGT");
        v.extend((0..49_996).map(|i| bases[i % bases.len()]));
        v
    };

    let pattern_medium = b"ACGTGCATGCTAGCTAGCTAGCTAGCTAGCATGCATGCATGCATGCATGC".as_slice();
    let pattern_long = b"ACGTGCATGCTAGCTAGCTAGCTAGCTAGCATGCATGCATGCATGCATGCACGTGCATGCTAGCTAGCTAGCTAGCTAGCATGCATGCATGCATGCATGCACGTGCATGCTAGCTAGCTAGCTAGCTAGC".as_slice();
    let pattern_very_long: Vec<u8> = {
        let bases = b"ACGTGCATGCTAGC";
        (0..500).map(|i| bases[i % bases.len()]).collect()
    };
    let long_medium_late: Vec<u8> = {
        let mut v = b"A".repeat(9_950);
        v.extend_from_slice(pattern_medium);
        v
    };
    let vlong_medium_late: Vec<u8> = {
        let mut v = b"A".repeat(99_950);
        v.extend_from_slice(pattern_medium);
        v
    };
    let vlong_long_late: Vec<u8> = {
        let mut v = b"A".repeat(99_870);
        v.extend_from_slice(pattern_long);
        v
    };
    let vlong_vlong_late: Vec<u8> = {
        let mut v = b"A".repeat(99_500);
        v.extend_from_slice(&pattern_very_long);
        v
    };

    // contains (4bp motif)
    {
        let filter = LikeFilter {
            op: SequenceOp::Contains,
            pattern: "ACGT".to_string(),
        };
        let cases: &[(&str, &[u8])] = &[
            ("short/no_match", &short_no_match),
            ("short/early_match", &short_early_match),
            ("short/late_match", &short_late_match),
            ("long/no_match", &long_no_match),
            ("long/early_match", &long_early_match),
            ("long/late_match", &long_late_match),
            ("long/multiple_matches", &long_multiple_matches),
            ("long/realistic", &long_realistic),
            ("long/near_miss", &long_near_miss),
            ("vlong/no_match", &vlong_no_match),
            ("vlong/late_match", &vlong_late_match),
            ("vlong/near_miss", &vlong_near_miss),
            ("vlong/realistic", &vlong_realistic),
        ];
        let mut group = c.benchmark_group("sequence_search/contains_4bp");
        for (name, seq) in cases {
            group.bench_function(format!("{name}"), |b| {
                b.iter(|| black_box(filter.eval(black_box(seq))))
            });
        }
        group.finish();
    }

    // contains (50bp probe)
    {
        let filter = LikeFilter {
            op: SequenceOp::Contains,
            pattern: String::from_utf8(pattern_medium.to_vec())
                .unwrap()
                .to_ascii_uppercase(),
        };
        let cases: &[(&str, &[u8])] = &[
            ("long/no_match", &long_no_match),
            ("long/late_match", &long_medium_late),
            ("long/near_miss", &long_near_miss),
            ("long/realistic", &long_realistic),
            ("vlong/no_match", &vlong_no_match),
            ("vlong/late_match", &vlong_medium_late),
            ("vlong/near_miss", &vlong_near_miss),
        ];
        let mut group = c.benchmark_group("sequence_search/contains_50bp");
        for (name, seq) in cases {
            group.bench_function(format!("{name}"), |b| {
                b.iter(|| black_box(filter.eval(black_box(seq))))
            });
        }
        group.finish();
    }

    // contains (130bp read-length)
    {
        let filter = LikeFilter {
            op: SequenceOp::Contains,
            pattern: String::from_utf8(pattern_long.to_vec())
                .unwrap()
                .to_ascii_uppercase(),
        };
        let cases: &[(&str, &[u8])] = &[
            ("long/no_match", &long_no_match),
            ("long/near_miss", &long_near_miss),
            ("vlong/no_match", &vlong_no_match),
            ("vlong/late_match", &vlong_long_late),
            ("vlong/near_miss", &vlong_near_miss),
        ];
        let mut group = c.benchmark_group("sequence_search/contains_130bp");
        for (name, seq) in cases {
            group.bench_function(format!("{name}"), |b| {
                b.iter(|| black_box(filter.eval(black_box(seq))))
            });
        }
        group.finish();
    }

    // contains (500bp gene-scale)
    {
        let filter = LikeFilter {
            op: SequenceOp::Contains,
            pattern: String::from_utf8(pattern_very_long.clone())
                .unwrap()
                .to_ascii_uppercase(),
        };
        let cases: &[(&str, &[u8])] = &[
            ("vlong/no_match", &vlong_no_match),
            ("vlong/late_match", &vlong_vlong_late),
            ("vlong/near_miss", &vlong_near_miss),
        ];
        let mut group = c.benchmark_group("sequence_search/contains_500bp");
        for (name, seq) in cases {
            group.bench_function(format!("{name}"), |b| {
                b.iter(|| black_box(filter.eval(black_box(seq))))
            });
        }
        group.finish();
    }

    // starts_with / ends_with / eq — O(m) already, included for completeness
    {
        let sw_filter = LikeFilter {
            op: SequenceOp::StartsWith,
            pattern: "ACGT".to_string(),
        };
        let ew_filter = LikeFilter {
            op: SequenceOp::EndsWith,
            pattern: "ACGT".to_string(),
        };
        let eq_filter = LikeFilter {
            op: SequenceOp::Eq,
            pattern: "ACGT".to_string(),
        };
        let cases: &[(&str, &[u8])] = &[
            ("short", &short_no_match),
            ("long", &long_no_match),
            ("vlong", &vlong_no_match),
        ];
        let mut group = c.benchmark_group("sequence_search/starts_ends_eq");
        for (name, seq) in cases {
            group.bench_function(format!("starts_with/{name}"), |b| {
                b.iter(|| black_box(sw_filter.eval(black_box(seq))))
            });
            group.bench_function(format!("ends_with/{name}"), |b| {
                b.iter(|| black_box(ew_filter.eval(black_box(seq))))
            });
            group.bench_function(format!("eq/{name}"), |b| {
                b.iter(|| black_box(eq_filter.eval(black_box(seq))))
            });
        }
        group.finish();
    }
}

fn bench_scalar_functions(c: &mut Criterion) {
    let short: Vec<u8> = b"ACGTACGTACGTACGT".to_vec();
    let long: Vec<u8> = {
        let bases = b"ACGT";
        (0..10_000).map(|i| bases[i % 4]).collect()
    };
    let long_lowercase: Vec<u8> = {
        let bases = b"acgt";
        (0..10_000).map(|i| bases[i % 4]).collect()
    };
    let long_with_n: Vec<u8> = (0..10_000)
        .map(|i| if i % 10 == 0 { b'N' } else { b'A' })
        .collect();
    let vlong: Vec<u8> = {
        let bases = b"ACGT";
        (0..100_000).map(|i| bases[i % 4]).collect()
    };

    // compute_gc
    let mut group = c.benchmark_group("compute_gc");
    for (name, seq) in [
        ("short", short.as_slice()),
        ("long", long.as_slice()),
        ("vlong", vlong.as_slice()),
    ] {
        group.bench_function(format!("{name}"), |b| {
            b.iter(|| black_box(sqlite_fastx::functions::compute_gc(black_box(seq))))
        });
    }
    group.finish();

    // n_count
    let mut group = c.benchmark_group("n_count");
    for (name, seq) in [
        ("short", short.as_slice()),
        ("long_with_n", long_with_n.as_slice()),
        ("vlong", vlong.as_slice()),
    ] {
        group.bench_function(format!("{name}"), |b| {
            b.iter(|| black_box(sqlite_fastx::functions::n_count(black_box(seq))))
        });
    }
    group.finish();

    // base_count
    let mut group = c.benchmark_group("base_count");
    for (name, seq) in [
        ("short", short.as_slice()),
        ("long", long.as_slice()),
        ("long_lowercase", long_lowercase.as_slice()),
        ("vlong", vlong.as_slice()),
    ] {
        group.bench_function(format!("{name}"), |b| {
            b.iter(|| black_box(sqlite_fastx::functions::base_count(black_box(seq), b'G')))
        });
    }
    group.finish();
}

// --- FAI vs. no-FAI benchmarks ---
// Run the generators first:
//   cargo run --bin gen_bench_wide
//   cargo run --bin gen_bench_reads

fn setup_fasta(path: &str) -> Database {
    let db = Database::open(":memory:").unwrap();
    init(&db).unwrap();
    db.execute(
        &format!("CREATE VIRTUAL TABLE fa USING fasta('{path}')"),
        (),
    )
    .unwrap();
    db
}

fn bench_wide_full_scan(c: &mut Criterion) {
    let db_fai = setup_fasta(BENCH_WIDE_FAI);
    let db_nofai = setup_fasta(BENCH_WIDE_NOFAI);
    let mut group = c.benchmark_group("wide/full_scan");
    group.bench_function("with_fai", |b| {
        b.iter(|| black_box(count(&db_fai, "SELECT COUNT(*) FROM fa")))
    });
    group.bench_function("no_fai", |b| {
        b.iter(|| black_box(count(&db_nofai, "SELECT COUNT(*) FROM fa")))
    });
    group.finish();
}

fn bench_wide_id_lookup(c: &mut Criterion) {
    let db_fai = setup_fasta(BENCH_WIDE_FAI);
    let db_nofai = setup_fasta(BENCH_WIDE_NOFAI);
    let sql_middle = format!("SELECT COUNT(*) FROM fa WHERE id = 'seq{}'", WIDE_COUNT / 2);
    let sql_last = format!("SELECT COUNT(*) FROM fa WHERE id = 'seq{}'", WIDE_COUNT - 1);
    let mut group = c.benchmark_group("wide/id_lookup");
    group.bench_function("middle/with_fai", |b| {
        b.iter(|| black_box(count(&db_fai, &sql_middle)))
    });
    group.bench_function("middle/no_fai", |b| {
        b.iter(|| black_box(count(&db_nofai, &sql_middle)))
    });
    group.bench_function("last/with_fai", |b| {
        b.iter(|| black_box(count(&db_fai, &sql_last)))
    });
    group.bench_function("last/no_fai", |b| {
        b.iter(|| black_box(count(&db_nofai, &sql_last)))
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_full_scan,
    bench_length_pushdown,
    bench_length_range,
    bench_sequence_contains,
    bench_combined,
    bench_sequence_search,
    bench_scalar_functions,
    bench_wide_full_scan,
    bench_wide_id_lookup,
);
criterion_main!(benches);
