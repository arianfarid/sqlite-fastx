use nucleodb::init;
use sqlite3_ext::{Database, FallibleIteratorMut, FromValue};

const TEST_FA: &str = "tests/fixtures/test.fa";
const TEST_FASTQ: &str = "tests/fixtures/test.fastq";

// Fixture records and their expected values:
//
// test.fa (6 records):
//   seq1  ACGT        len=4   gc=0.50
//   seq2  AAAA        len=4   gc=0.00
//   seq3  GCGCGC      len=6   gc=1.00
//   seq4  ACGTTTTT    len=8   gc=0.25
//   seq5  TTTTACGT    len=8   gc=0.25
//   seq6  TTTTTTTTTT  len=10  gc=0.00
//
// test.fastq (3 records):
//   read1  ACGT      len=4  gc=0.50  quality=IIII (Q40)
//   read2  GGGG      len=4  gc=1.00  quality=!!!! (Q0)
//   read3  ACGTTTTT  len=8  gc=0.25  quality=???????? (Q30)

fn db() -> Database {
    let db = Database::open(":memory:").unwrap();
    init(&db).unwrap();
    db
}

fn fasta_db() -> Database {
    let db = db();
    db.execute(
        &format!("CREATE VIRTUAL TABLE fa USING fasta('{TEST_FA}')"),
        (),
    )
    .unwrap();
    db
}

fn fastq_db() -> Database {
    let db = db();
    db.execute(
        &format!("CREATE VIRTUAL TABLE fq USING fastq('{TEST_FASTQ}')"),
        (),
    )
    .unwrap();
    db
}

fn scalar_i64(db: &Database, sql: &str) -> i64 {
    let mut stmt = db.prepare(sql).unwrap();
    let rows = stmt.query(()).unwrap();
    let row = rows.next().unwrap().unwrap();
    row[0].get_i64()
}

fn scalar_f64(db: &Database, sql: &str) -> f64 {
    let mut stmt = db.prepare(sql).unwrap();
    let rows = stmt.query(()).unwrap();
    let row = rows.next().unwrap().unwrap();
    row[0].get_f64()
}

fn scalar_str(db: &Database, sql: &str) -> String {
    let mut stmt = db.prepare(sql).unwrap();
    let rows = stmt.query(()).unwrap();
    let row = rows.next().unwrap().unwrap();
    row[0].get_str().unwrap().to_string()
}

fn try_query(db: &Database, sql: &str) -> sqlite3_ext::Result<()> {
    let mut stmt = db.prepare(sql)?;
    let rows = stmt.query(())?;
    rows.next()?;
    Ok(())
}

// --- gc_content ---

#[test]
fn gc_content_half() {
    assert_eq!(scalar_f64(&db(), "SELECT gc_content('ACGT')"), 0.5);
}

#[test]
fn gc_content_empty() {
    assert_eq!(scalar_f64(&db(), "SELECT gc_content('')"), 0.0);
}

#[test]
fn gc_content_all_gc() {
    assert_eq!(scalar_f64(&db(), "SELECT gc_content('GCGC')"), 1.0);
}

#[test]
fn gc_content_none() {
    assert_eq!(scalar_f64(&db(), "SELECT gc_content('AAAA')"), 0.0);
}

// --- n_count ---

#[test]
fn n_count_basic() {
    assert_eq!(scalar_i64(&db(), "SELECT n_count('ACNGNTN')"), 3);
}

#[test]
fn n_count_none() {
    assert_eq!(scalar_i64(&db(), "SELECT n_count('ACGT')"), 0);
}

#[test]
fn n_count_lowercase() {
    assert_eq!(scalar_i64(&db(), "SELECT n_count('acngntn')"), 3);
}

// --- base_count ---

#[test]
fn base_count_single() {
    assert_eq!(scalar_i64(&db(), "SELECT base_count('ACGT', 'A')"), 1);
}

#[test]
fn base_count_repeated() {
    assert_eq!(scalar_i64(&db(), "SELECT base_count('AAAA', 'A')"), 4);
}

#[test]
fn base_count_case_insensitive() {
    assert_eq!(scalar_i64(&db(), "SELECT base_count('AaCcGgTt', 'a')"), 2);
}

#[test]
fn base_count_invalid_base_errors() {
    assert!(try_query(&db(), "SELECT base_count('ACGT', 'Z')").is_err());
}

// --- to_rna / to_dna ---

#[test]
fn to_rna_converts_t() {
    assert_eq!(scalar_str(&db(), "SELECT to_rna('ACGT')"), "ACGU");
}

#[test]
fn to_rna_lowercase() {
    assert_eq!(scalar_str(&db(), "SELECT to_rna('acgt')"), "acgu");
}

#[test]
fn to_dna_converts_u() {
    assert_eq!(scalar_str(&db(), "SELECT to_dna('ACGU')"), "ACGT");
}

#[test]
fn to_dna_lowercase() {
    assert_eq!(scalar_str(&db(), "SELECT to_dna('acgu')"), "acgt");
}

// --- reverse_complement ---

#[test]
fn reverse_complement_palindrome() {
    assert_eq!(scalar_str(&db(), "SELECT reverse_complement('ACGT')"), "ACGT");
}

#[test]
fn reverse_complement_asymmetric() {
    assert_eq!(
        scalar_str(&db(), "SELECT reverse_complement('AAAA')"),
        "TTTT"
    );
}

#[test]
fn reverse_complement_roundtrip() {
    assert_eq!(
        scalar_str(
            &db(),
            "SELECT reverse_complement(reverse_complement('ACGTACGT'))"
        ),
        "ACGTACGT"
    );
}

// --- is_valid_dna / is_valid_rna ---

#[test]
fn is_valid_dna_accepts_acgtn() {
    assert_eq!(scalar_i64(&db(), "SELECT is_valid_dna('ACGTN')"), 1);
}

#[test]
fn is_valid_dna_rejects_u() {
    assert_eq!(scalar_i64(&db(), "SELECT is_valid_dna('ACGTU')"), 0);
}

#[test]
fn is_valid_dna_rejects_invalid() {
    assert_eq!(scalar_i64(&db(), "SELECT is_valid_dna('ACGTZ')"), 0);
}

#[test]
fn is_valid_rna_accepts_acgun() {
    assert_eq!(scalar_i64(&db(), "SELECT is_valid_rna('ACGUN')"), 1);
}

#[test]
fn is_valid_rna_rejects_t() {
    assert_eq!(scalar_i64(&db(), "SELECT is_valid_rna('ACGUT')"), 0);
}

// --- min_quality / mean_quality ---

#[test]
fn min_quality_uniform_high() {
    // 'I' = ASCII 73, Phred 40
    assert_eq!(scalar_i64(&db(), "SELECT min_quality('IIII')"), 40);
}

#[test]
fn min_quality_finds_minimum() {
    // '!' = Q0 is the minimum among the I's
    assert_eq!(scalar_i64(&db(), "SELECT min_quality('III!')"), 0);
}

#[test]
fn mean_quality_uniform() {
    assert_eq!(scalar_f64(&db(), "SELECT mean_quality('IIII')"), 40.0);
}

#[test]
fn mean_quality_mixed() {
    // '!' = Q0, 'I' = Q40, mean = 20.0
    assert_eq!(scalar_f64(&db(), "SELECT mean_quality('!I')"), 20.0);
}

// --- n50 ---

#[test]
fn n50_five_contigs() {
    // total=1500, half=750; sorted desc: 500,400,300,200,100
    // running sum: 500 (< 750), 900 (>= 750) -> N50 = 400
    assert_eq!(
        scalar_i64(
            &db(),
            "SELECT n50(column1) FROM (VALUES (100), (200), (300), (400), (500))"
        ),
        400
    );
}

#[test]
fn n50_single_contig() {
    assert_eq!(
        scalar_i64(&db(), "SELECT n50(column1) FROM (VALUES (42))"),
        42
    );
}

#[test]
fn n50_on_fasta_table() {
    // lengths: 4,4,6,8,8,10; total=40, half=20
    // sorted desc: 10,8,8,6,4,4; running: 10->18->26 >= 20 -> N50 = 8
    assert_eq!(scalar_i64(&fasta_db(), "SELECT n50(length) FROM fa"), 8);
}

// --- FASTA virtual table ---

#[test]
fn fasta_row_count() {
    assert_eq!(scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa"), 6);
}

#[test]
fn fasta_columns() {
    let db = fasta_db();
    assert_eq!(scalar_str(&db, "SELECT id FROM fa WHERE id = 'seq1'"), "seq1");
    assert_eq!(
        scalar_str(&db, "SELECT sequence FROM fa WHERE id = 'seq1'"),
        "ACGT"
    );
    assert_eq!(scalar_i64(&db, "SELECT length FROM fa WHERE id = 'seq1'"), 4);
    assert_eq!(
        scalar_f64(&db, "SELECT gc_content FROM fa WHERE id = 'seq1'"),
        0.5
    );
}

#[test]
fn fasta_length_gt() {
    // seq4(8), seq5(8), seq6(10)
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE length > 6"),
        3
    );
}

#[test]
fn fasta_length_ge() {
    // seq3(6), seq4(8), seq5(8), seq6(10)
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE length >= 6"),
        4
    );
}

#[test]
fn fasta_length_lt() {
    // seq1(4), seq2(4)
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE length < 6"),
        2
    );
}

#[test]
fn fasta_length_eq() {
    // seq1(4), seq2(4)
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE length = 4"),
        2
    );
}

#[test]
fn fasta_length_range() {
    // seq3(6), seq4(8), seq5(8)
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE length > 4 AND length < 10"
        ),
        3
    );
}

#[test]
fn fasta_sequence_contains() {
    // seq1(ACGT), seq4(ACGTTTTT), seq5(TTTTACGT) all contain the motif
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE sequence LIKE '%ACGT%'"
        ),
        3
    );
}

#[test]
fn fasta_sequence_starts_with() {
    // seq4(ACGTTTTT) and seq1(ACGT) start with ACGT; seq5(TTTTACGT) does not
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE sequence LIKE 'ACGT%'"
        ),
        2
    );
}

#[test]
fn fasta_sequence_ends_with() {
    // seq5(TTTTACGT) and seq1(ACGT) end with ACGT; seq4(ACGTTTTT) does not
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE sequence LIKE '%ACGT'"
        ),
        2
    );
}

#[test]
fn fasta_sequence_exact() {
    // only seq1 is exactly ACGT
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE sequence LIKE 'ACGT'"
        ),
        1
    );
}

#[test]
fn fasta_gc_content_gt() {
    // seq1(0.50), seq3(1.00) are above 0.4; seq4/seq5 are 0.25, seq2/seq6 are 0.0
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE gc_content > 0.4"
        ),
        2
    );
}

#[test]
fn fasta_file_not_found() {
    let db = db();
    db.execute(
        "CREATE VIRTUAL TABLE bad USING fasta('nonexistent.fa')",
        (),
    )
    .unwrap();
    assert!(try_query(&db, "SELECT * FROM bad").is_err());
}

// --- FASTQ virtual table ---

#[test]
fn fastq_row_count() {
    assert_eq!(scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq"), 3);
}

#[test]
fn fastq_quality_column() {
    assert_eq!(
        scalar_str(&fastq_db(), "SELECT quality FROM fq WHERE id = 'read1'"),
        "IIII"
    );
}

#[test]
fn fastq_length_filter() {
    // only read3 has length 8 (> 4)
    assert_eq!(
        scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq WHERE length > 4"),
        1
    );
}

#[test]
fn fastq_sequence_contains() {
    // read1(ACGT) and read3(ACGTTTTT) both contain ACGT
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE sequence LIKE '%ACGT%'"
        ),
        2
    );
}

#[test]
fn fastq_file_not_found() {
    let db = db();
    db.execute(
        "CREATE VIRTUAL TABLE bad USING fastq('nonexistent.fastq')",
        (),
    )
    .unwrap();
    assert!(try_query(&db, "SELECT * FROM bad").is_err());
}
