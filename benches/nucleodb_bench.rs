use criterion::{Criterion, criterion_group, criterion_main};
use nucleodb::init;
use sqlite3_ext::{Database, FallibleIteratorMut, FromValue};
use std::hint::black_box;

const BENCH_FA: &str = "bench.fa";

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
        b.iter(|| black_box(count(&db, "SELECT COUNT(*) FROM fa WHERE length + 0 > 100 AND CAST(sequence AS TEXT) LIKE '%ACGT%'")))
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
);
criterion_main!(benches);
