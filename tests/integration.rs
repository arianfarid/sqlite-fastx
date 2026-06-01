use sqlite_fastx::init;
use sqlite3_ext::{Database, FallibleIteratorMut, FromValue};

const TEST_FA: &str = "tests/fixtures/test.fa";
const TEST_FASTQ: &str = "tests/fixtures/test.fastq";
const TEST_FASTQ_FAI: &str = "tests/fixtures/test.fastq.fai";
const TEST_BGZF_FA: &str = "tests/fixtures/sample.fa.gz";
const TEST_BGZF_FAI: &str = "tests/fixtures/sample.fa.gz.fai";
const TEST_BGZF_FASTQ: &str = "tests/fixtures/sample.fastq.gz";
const TEST_BGZF_FASTQ_FAI: &str = "tests/fixtures/sample.fastq.gz.fai";

// Fixture records and their expected values:
//
// test.fa (8 records):
//   seq1  ACGT        len=4   gc=0.50  desc="four bases half gc"
//   seq2  AAAA        len=4   gc=0.00  desc="four adenine no gc"
//   seq3  GCGCGC      len=6   gc=1.00  desc="six bases pure gc"
//   seq4  ACGTTTTT    len=8   gc=0.25  desc="eight bases starts with acgt"
//   seq5  TTTTACGT    len=8   gc=0.25  desc="eight bases ends with acgt"
//   seq6  TTTTTTTTTT  len=10  gc=0.00  desc="ten thymine no gc"
//   seq7  acgtacgt    len=8   gc=0.50  desc="lowercase sequence mixed case"
//   seq8  ACGNNNTT    len=8   gc=0.25  desc="ambiguous bases with n"
//
// test.fastq (4 records):
//   read1  ACGT      len=4  gc=0.50  quality=IIII     (Q40)  desc="high quality"
//   read2  GGGG      len=4  gc=1.00  quality=!!!!     (Q0)   desc="low quality"
//   read3  ACGTTTTT  len=8  gc=0.25  quality=???????? (Q30)  desc="medium quality"
//   read4  acgtgggg  len=8  gc=0.75  quality=IIII???? (mixed) desc="lowercase sequence"

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

// --- gc_content scalar ---

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

// --- n_count scalar ---

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

// --- base_count scalar ---

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

// --- to_rna / to_dna scalar ---

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

// --- to_rna / to_dna: no-op cases ---

#[test]
fn to_rna_already_rna() {
    assert_eq!(scalar_str(&db(), "SELECT to_rna('ACGU')"), "ACGU");
}

#[test]
fn to_dna_already_dna() {
    assert_eq!(scalar_str(&db(), "SELECT to_dna('ACGT')"), "ACGT");
}

// --- reverse_complement scalar ---

#[test]
fn reverse_complement_palindrome() {
    assert_eq!(
        scalar_str(&db(), "SELECT reverse_complement('ACGT')"),
        "ACGT"
    );
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

// --- complement scalar ---

#[test]
fn complement_single_bases() {
    assert_eq!(scalar_str(&db(), "SELECT complement('A')"), "T");
    assert_eq!(scalar_str(&db(), "SELECT complement('T')"), "A");
    assert_eq!(scalar_str(&db(), "SELECT complement('G')"), "C");
    assert_eq!(scalar_str(&db(), "SELECT complement('C')"), "G");
    assert_eq!(scalar_str(&db(), "SELECT complement('U')"), "A");
}

#[test]
fn complement_empty() {
    assert_eq!(scalar_str(&db(), "SELECT complement('')"), "");
}

#[test]
fn complement_does_not_reverse() {
    assert_eq!(scalar_str(&db(), "SELECT complement('AAAC')"), "TTTG");
}

#[test]
fn complement_lowercase() {
    assert_eq!(scalar_str(&db(), "SELECT complement('acgt')"), "tgca");
}

#[test]
fn complement_mixed_case() {
    assert_eq!(scalar_str(&db(), "SELECT complement('AcGt')"), "TgCa");
}

#[test]
fn complement_ambiguous_passthrough() {
    assert_eq!(scalar_str(&db(), "SELECT complement('N')"), "N");
    assert_eq!(scalar_str(&db(), "SELECT complement('n')"), "n");
}

#[test]
fn complement_rna() {
    assert_eq!(scalar_str(&db(), "SELECT complement('ACGU')"), "TGCA");
}

#[test]
fn complement_differs_from_reverse_complement() {
    // complement('AAAC') = 'TTTG', reverse_complement('AAAC') = 'GTTT'
    assert_ne!(
        scalar_str(&db(), "SELECT complement('AAAC')"),
        scalar_str(&db(), "SELECT reverse_complement('AAAC')")
    );
}

#[test]
fn complement_rc_roundtrip() {
    // complement then reverse equals reverse_complement
    assert_eq!(
        scalar_str(
            &db(),
            "SELECT reverse_complement('ACGTACGT')"
        ),
        scalar_str(
            &db(),
            "SELECT complement('ACGTACGT')"
        ).chars().rev().collect::<String>()
    );
}

// --- is_valid_dna / is_valid_rna scalar ---

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

// --- min_quality / mean_quality scalar ---

#[test]
fn min_quality_uniform_high() {
    assert_eq!(scalar_i64(&db(), "SELECT min_quality('IIII')"), 40);
}

#[test]
fn min_quality_finds_minimum() {
    assert_eq!(scalar_i64(&db(), "SELECT min_quality('III!')"), 0);
}

#[test]
fn mean_quality_uniform() {
    assert_eq!(scalar_f64(&db(), "SELECT mean_quality('IIII')"), 40.0);
}

#[test]
fn mean_quality_mixed() {
    assert_eq!(scalar_f64(&db(), "SELECT mean_quality('!I')"), 20.0);
}

// --- max_quality scalar ---

#[test]
fn max_quality_uniform_high() {
    assert_eq!(scalar_i64(&db(), "SELECT max_quality('IIII')"), 40);
}

#[test]
fn max_quality_finds_maximum() {
    assert_eq!(scalar_i64(&db(), "SELECT max_quality('!!!I')"), 40);
}

#[test]
fn max_quality_all_low() {
    assert_eq!(scalar_i64(&db(), "SELECT max_quality('!!!!')"), 0);
}

#[test]
fn max_quality_empty() {
    assert_eq!(scalar_i64(&db(), "SELECT max_quality('')"), 0);
}

// --- n50 aggregate ---

#[test]
fn n50_five_contigs() {
    // total=1500, half=750; sorted desc: 500,400,300,200,100
    // running: 500 < 750, 900 >= 750 -> N50 = 400
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
    // lengths: 4,4,6,8,8,8,8,10; total=56, half=28
    // sorted desc: 10,8,8,8,8,6,4,4
    // running: 10->18->26->34 >= 28 -> N50 = 8
    assert_eq!(scalar_i64(&fasta_db(), "SELECT n50(length) FROM fa"), 8);
}

// --- base_totals aggregate ---

#[test]
fn base_totals_pure_dna() {
    // "AACGT" -> A:2, C:1, G:1, T:1, U:0
    assert_eq!(
        scalar_str(
            &db(),
            "SELECT base_totals(column1) FROM (VALUES ('AACGT'))"
        ),
        r#"{"A": 2, "C": 1, "G": 1, "T": 1, "U": 0}"#
    );
}

#[test]
fn base_totals_accumulates_across_rows() {
    // "AAAA" + "CCCC" -> A:4, C:4, G:0, T:0, U:0
    assert_eq!(
        scalar_str(
            &db(),
            "SELECT base_totals(column1) FROM (VALUES ('AAAA'), ('CCCC'))"
        ),
        r#"{"A": 4, "C": 4, "G": 0, "T": 0, "U": 0}"#
    );
}

#[test]
fn base_totals_rna() {
    // U counted separately, T not present
    assert_eq!(
        scalar_str(
            &db(),
            "SELECT base_totals(column1) FROM (VALUES ('ACGU'))"
        ),
        r#"{"A": 1, "C": 1, "G": 1, "T": 0, "U": 1}"#
    );
}

#[test]
fn base_totals_n_excluded() {
    // N does not count toward any bucket
    assert_eq!(
        scalar_str(
            &db(),
            "SELECT base_totals(column1) FROM (VALUES ('ACGTN'))"
        ),
        r#"{"A": 1, "C": 1, "G": 1, "T": 1, "U": 0}"#
    );
}

#[test]
fn base_totals_lowercase() {
    assert_eq!(
        scalar_str(
            &db(),
            "SELECT base_totals(column1) FROM (VALUES ('acgt'))"
        ),
        r#"{"A": 1, "C": 1, "G": 1, "T": 1, "U": 0}"#
    );
}

#[test]
fn base_totals_empty_row() {
    assert_eq!(
        scalar_str(
            &db(),
            "SELECT base_totals(column1) FROM (VALUES (''))"
        ),
        r#"{"A": 0, "C": 0, "G": 0, "T": 0, "U": 0}"#
    );
}

#[test]
fn base_totals_on_fasta_table() {
    // test.fa has 8 records; just verify A count is a positive integer
    let result = scalar_str(&fasta_db(), "SELECT base_totals(sequence) FROM fa");
    assert!(result.contains("\"A\":"));
    let a_val: f64 = result
        .split("\"A\": ")
        .nth(1)
        .unwrap()
        .split(',')
        .next()
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert!(a_val > 0.0);
}

// --- FASTA virtual table: basic ---

#[test]
fn fasta_row_count() {
    assert_eq!(scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa"), 8);
}

#[test]
fn fasta_columns() {
    let db = fasta_db();
    assert_eq!(
        scalar_str(&db, "SELECT id FROM fa WHERE id = 'seq1'"),
        "seq1"
    );
    assert_eq!(
        scalar_str(&db, "SELECT sequence FROM fa WHERE id = 'seq1'"),
        "ACGT"
    );
    assert_eq!(
        scalar_i64(&db, "SELECT length FROM fa WHERE id = 'seq1'"),
        4
    );
    assert_eq!(
        scalar_f64(&db, "SELECT gc_content FROM fa WHERE id = 'seq1'"),
        0.5
    );
}

#[test]
fn fasta_id_column_nonexistent_returns_zero_rows() {
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE id = 'notreal'"),
        0
    );
}

#[test]
fn fasta_description_column() {
    assert_eq!(
        scalar_str(&fasta_db(), "SELECT description FROM fa WHERE id = 'seq1'"),
        "four bases half gc"
    );
}

// --- FASTA virtual table: length filters ---

#[test]
fn fasta_length_gt() {
    // seq4(8), seq5(8), seq6(10), seq7(8), seq8(8)
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE length > 6"),
        5
    );
}

#[test]
fn fasta_length_ge() {
    // seq3(6), seq4(8), seq5(8), seq6(10), seq7(8), seq8(8)
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE length >= 6"),
        6
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
fn fasta_length_le() {
    // seq1(4), seq2(4), seq3(6)
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE length <= 6"),
        3
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
    // seq3(6), seq4(8), seq5(8), seq7(8), seq8(8)
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE length > 4 AND length < 10"
        ),
        5
    );
}

// --- FASTA virtual table: sequence LIKE filters ---

#[test]
fn fasta_sequence_contains() {
    // seq1(ACGT), seq4(ACGTTTTT), seq5(TTTTACGT), seq7(acgtacgt lowercase)
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE sequence LIKE '%ACGT%'"
        ),
        4
    );
}

#[test]
fn fasta_sequence_starts_with() {
    // seq1(ACGT), seq4(ACGTTTTT), seq7(acgtacgt)
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE sequence LIKE 'ACGT%'"
        ),
        3
    );
}

#[test]
fn fasta_sequence_ends_with() {
    // seq1(ACGT), seq5(TTTTACGT), seq7(acgtacgt)
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE sequence LIKE '%ACGT'"
        ),
        3
    );
}

#[test]
fn fasta_sequence_exact() {
    // only seq1 is exactly ACGT (seq7 is acgtacgt, different length)
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE sequence LIKE 'ACGT'"
        ),
        1
    );
}

#[test]
fn fasta_sequence_lowercase_matches_uppercase_pattern() {
    // seq7 has lowercase sequence acgtacgt, should match uppercase pattern
    assert_eq!(
        scalar_str(
            &fasta_db(),
            "SELECT id FROM fa WHERE sequence LIKE 'acgtacgt'"
        ),
        "seq7"
    );
}

// --- FASTA virtual table: gc_content filters ---

#[test]
fn fasta_gc_content_gt() {
    // seq1(0.50), seq3(1.00), seq7(0.50) above 0.4
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE gc_content > 0.4"
        ),
        3
    );
}

#[test]
fn fasta_gc_content_ge() {
    // seq1(0.50), seq3(1.00), seq7(0.50) >= 0.5
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE gc_content >= 0.5"
        ),
        3
    );
}

#[test]
fn fasta_gc_content_lt() {
    // seq2(0.0), seq4(0.25), seq5(0.25), seq6(0.0), seq8(0.25) < 0.5
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE gc_content < 0.5"
        ),
        5
    );
}

#[test]
fn fasta_gc_content_le() {
    // everything except seq3(1.0)
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE gc_content <= 0.5"
        ),
        7
    );
}

// --- FASTA virtual table: id LIKE filters ---

#[test]
fn fasta_id_starts_with() {
    assert_eq!(
        scalar_i64(&fasta_db(), "SELECT COUNT(*) FROM fa WHERE id LIKE 'seq%'"),
        8
    );
}

#[test]
fn fasta_id_exact() {
    assert_eq!(
        scalar_str(&fasta_db(), "SELECT id FROM fa WHERE id LIKE 'seq1'"),
        "seq1"
    );
}

#[test]
fn fasta_id_case_insensitive() {
    assert_eq!(
        scalar_str(&fasta_db(), "SELECT id FROM fa WHERE id LIKE 'SEQ1'"),
        "seq1"
    );
}

#[test]
fn fasta_id_ends_with() {
    // seq1 through seq9 — only seq1 ends with '1'... wait, just seq1 and no others end with same digit
    // seq1 ends with '1', use that
    assert_eq!(
        scalar_str(&fasta_db(), "SELECT id FROM fa WHERE id LIKE '%1'"),
        "seq1"
    );
}

// --- FASTA virtual table: description LIKE filters ---

#[test]
fn fasta_description_contains() {
    // seq4 "eight bases starts with acgt", seq5 "eight bases ends with acgt"
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE description LIKE '%acgt%'"
        ),
        2
    );
}

#[test]
fn fasta_description_case_insensitive() {
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE description LIKE '%ACGT%'"
        ),
        2
    );
}

#[test]
fn fasta_description_starts_with() {
    // seq1, seq2 both start with "four"
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE description LIKE 'four%'"
        ),
        2
    );
}

// --- FASTA virtual table: combined filters ---

#[test]
fn fasta_length_and_gc_combined() {
    // length > 4 AND gc_content >= 0.5: seq3(6,1.0), seq7(8,0.5)
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE length > 4 AND gc_content >= 0.5"
        ),
        2
    );
}

#[test]
fn fasta_id_and_length_combined() {
    // id LIKE 'seq%' AND length = 4: seq1, seq2
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE id LIKE 'seq%' AND length = 4"
        ),
        2
    );
}

// --- FASTA virtual table: n_count on fixture ---

#[test]
fn fasta_n_count_on_record() {
    // seq8 has 3 N bases
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT n_count(sequence) FROM fa WHERE id = 'seq8'"
        ),
        3
    );
}

#[test]
fn fasta_n_count_zero_on_clean_record() {
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT n_count(sequence) FROM fa WHERE id = 'seq1'"
        ),
        0
    );
}

// --- FASTA virtual table: error handling ---

#[test]
fn fasta_file_not_found() {
    let db = db();
    db.execute("CREATE VIRTUAL TABLE bad USING fasta('nonexistent.fa')", ())
        .unwrap();
    assert!(try_query(&db, "SELECT * FROM bad").is_err());
}

// --- FASTQ virtual table: basic ---

#[test]
fn fastq_row_count() {
    assert_eq!(scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq"), 4);
}

#[test]
fn fastq_quality_column() {
    assert_eq!(
        scalar_str(&fastq_db(), "SELECT quality FROM fq WHERE id = 'read1'"),
        "IIII"
    );
}

#[test]
fn fastq_description_column() {
    assert_eq!(
        scalar_str(&fastq_db(), "SELECT description FROM fq WHERE id = 'read1'"),
        "high quality"
    );
}

// --- FASTQ virtual table: length filters ---

#[test]
fn fastq_length_filter() {
    // read3(8) and read4(8) have length > 4
    assert_eq!(
        scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq WHERE length > 4"),
        2
    );
}

// --- FASTQ virtual table: sequence LIKE filters ---

#[test]
fn fastq_sequence_contains() {
    // read1(ACGT), read3(ACGTTTTT), read4(acgtgggg lowercase)
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE sequence LIKE '%ACGT%'"
        ),
        3
    );
}

#[test]
fn fastq_sequence_lowercase_matches() {
    assert_eq!(
        scalar_str(
            &fastq_db(),
            "SELECT id FROM fq WHERE sequence LIKE 'acgtgggg'"
        ),
        "read4"
    );
}

// --- FASTQ virtual table: gc_content filters ---

#[test]
fn fastq_gc_content_gt() {
    // read1(0.50), read2(1.00), read4(0.75) > 0.4; read3(0.25) not
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE gc_content > 0.4"
        ),
        3
    );
}

#[test]
fn fastq_gc_content_ge() {
    // read1(0.50), read2(1.00), read4(0.75) >= 0.5
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE gc_content >= 0.5"
        ),
        3
    );
}

// --- FASTQ virtual table: id LIKE filters ---

#[test]
fn fastq_id_starts_with() {
    assert_eq!(
        scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq WHERE id LIKE 'read%'"),
        4
    );
}

#[test]
fn fastq_id_exact_case_insensitive() {
    assert_eq!(
        scalar_str(&fastq_db(), "SELECT id FROM fq WHERE id LIKE 'READ1'"),
        "read1"
    );
}

// --- FASTQ virtual table: description LIKE filters ---

#[test]
fn fastq_description_contains() {
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE description LIKE '%quality%'"
        ),
        3 // high quality, low quality, medium quality
    );
}

#[test]
fn fastq_description_exact() {
    assert_eq!(
        scalar_str(
            &fastq_db(),
            "SELECT id FROM fq WHERE description LIKE 'lowercase sequence'"
        ),
        "read4"
    );
}

// --- FASTQ virtual table: error handling ---

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

// --- FASTQ virtual table: missing length filter variants ---

#[test]
fn fastq_length_ge() {
    // all four reads have length >= 4
    assert_eq!(
        scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq WHERE length >= 4"),
        4
    );
}

#[test]
fn fastq_length_lt() {
    // read1(4), read2(4) < 8
    assert_eq!(
        scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq WHERE length < 8"),
        2
    );
}

#[test]
fn fastq_length_le() {
    // all four reads have length <= 8
    assert_eq!(
        scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq WHERE length <= 8"),
        4
    );
}

#[test]
fn fastq_length_eq() {
    // read1(4), read2(4)
    assert_eq!(
        scalar_i64(&fastq_db(), "SELECT COUNT(*) FROM fq WHERE length = 4"),
        2
    );
}

// --- FASTQ virtual table: gc_content lt/le ---

#[test]
fn fastq_gc_content_lt() {
    // read3(0.25) < 0.5
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE gc_content < 0.5"
        ),
        1
    );
}

#[test]
fn fastq_gc_content_le() {
    // read1(0.50), read3(0.25) <= 0.5
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE gc_content <= 0.5"
        ),
        2
    );
}

// --- FASTQ virtual table: mean_quality / min_quality column values ---

#[test]
fn fastq_mean_quality_column_high() {
    // read1: IIII → Q40 each → mean = 40.0
    assert_eq!(
        scalar_f64(
            &fastq_db(),
            "SELECT mean_quality FROM fq WHERE id LIKE 'read1'"
        ),
        40.0
    );
}

#[test]
fn fastq_mean_quality_column_mixed() {
    // read4: IIII???? → Q40,Q40,Q40,Q40,Q30,Q30,Q30,Q30 → mean = 35.0
    assert_eq!(
        scalar_f64(
            &fastq_db(),
            "SELECT mean_quality FROM fq WHERE id LIKE 'read4'"
        ),
        35.0
    );
}

#[test]
fn fastq_min_quality_column_high() {
    // read1: IIII → all Q40 → min = 40
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT min_quality FROM fq WHERE id LIKE 'read1'"
        ),
        40
    );
}

#[test]
fn fastq_min_quality_column_mixed() {
    // read4: IIII???? → Q40 and Q30 mixed → min = 30
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT min_quality FROM fq WHERE id LIKE 'read4'"
        ),
        30
    );
}

// --- FASTQ virtual table: mean_quality filter ---

#[test]
fn fastq_mean_quality_gt() {
    // read1(40.0), read4(35.0) > 30
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE mean_quality > 30"
        ),
        2
    );
}

#[test]
fn fastq_mean_quality_ge() {
    // read1(40.0), read3(30.0), read4(35.0) >= 30
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE mean_quality >= 30"
        ),
        3
    );
}

#[test]
fn fastq_mean_quality_lt() {
    // read2(0.0) < 30
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE mean_quality < 30"
        ),
        1
    );
}

#[test]
fn fastq_mean_quality_le() {
    // read2(0.0), read3(30.0) <= 30
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE mean_quality <= 30"
        ),
        2
    );
}

// --- FASTQ virtual table: min_quality filter ---

#[test]
fn fastq_min_quality_gt() {
    // read1(40) > 30
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE min_quality > 30"
        ),
        1
    );
}

#[test]
fn fastq_min_quality_ge() {
    // read1(40), read3(30), read4(30) >= 30
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE min_quality >= 30"
        ),
        3
    );
}

#[test]
fn fastq_min_quality_lt() {
    // read2(0) < 30
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE min_quality < 30"
        ),
        1
    );
}

// --- FASTQ virtual table: combined filters ---

#[test]
fn fastq_length_and_mean_quality_combined() {
    // length > 4 AND mean_quality >= 30: read3(8, 30.0), read4(8, 35.0)
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE length > 4 AND mean_quality >= 30"
        ),
        2
    );
}

#[test]
fn fastq_gc_content_and_min_quality_combined() {
    // gc_content > 0.5 AND min_quality > 0:
    //   read2(gc=1.0, min=0) → excluded (min not > 0)
    //   read4(gc=0.75, min=30) → included
    assert_eq!(
        scalar_i64(
            &fastq_db(),
            "SELECT COUNT(*) FROM fq WHERE gc_content > 0.5 AND min_quality > 0"
        ),
        1
    );
}

// --- FASTA virtual table: id equality (Eq path, sets unique=true) ---

#[test]
fn fasta_id_eq() {
    assert_eq!(
        scalar_str(&fasta_db(), "SELECT id FROM fa WHERE id = 'seq1'"),
        "seq1"
    );
}

// --- base_composition scalar ---

#[test]
fn base_composition_equal_dna() {
    assert_eq!(
        scalar_str(&db(), "SELECT base_composition('ACGT')"),
        r#"{"A": 0.2500, "C": 0.2500, "G": 0.2500, "T": 0.2500, "U": 0.0000}"#
    );
}

#[test]
fn base_composition_empty() {
    assert_eq!(
        scalar_str(&db(), "SELECT base_composition('')"),
        r#"{"A": 0.0000, "C": 0.0000, "G": 0.0000, "T": 0.0000, "U": 0.0000}"#
    );
}

#[test]
fn base_composition_pure_gc() {
    assert_eq!(
        scalar_str(&db(), "SELECT base_composition('GCGC')"),
        r#"{"A": 0.0000, "C": 0.5000, "G": 0.5000, "T": 0.0000, "U": 0.0000}"#
    );
}

#[test]
fn base_composition_rna() {
    assert_eq!(
        scalar_str(&db(), "SELECT base_composition('ACGU')"),
        r#"{"A": 0.2500, "C": 0.2500, "G": 0.2500, "T": 0.0000, "U": 0.2500}"#
    );
}

// --- gc_content scalar: ambiguous bases ---

#[test]
fn gc_content_with_n_bases() {
    // ACGTN: 2 GC out of 5 = 0.4
    assert_eq!(scalar_f64(&db(), "SELECT gc_content('ACGTN')"), 0.4);
}

// --- reverse_complement: lowercase via SQL ---

#[test]
fn reverse_complement_lowercase_sql() {
    assert_eq!(
        scalar_str(&db(), "SELECT reverse_complement('aaaa')"),
        "tttt"
    );
}

// --- is_valid_dna / is_valid_rna: empty string ---

#[test]
fn is_valid_dna_empty() {
    assert_eq!(scalar_i64(&db(), "SELECT is_valid_dna('')"), 1);
}

#[test]
fn is_valid_rna_empty() {
    assert_eq!(scalar_i64(&db(), "SELECT is_valid_rna('')"), 1);
}

// --- fastq fai seek ---

fn fq_fai_db() -> Database {
    assert!(
        std::path::Path::new(TEST_FASTQ_FAI).exists(),
        "FAI fixture missing — run: samtools faidx {TEST_FASTQ}"
    );
    let db = db();
    db.execute(
        &format!("CREATE VIRTUAL TABLE fq_fai USING fastq('{TEST_FASTQ}')"),
        (),
    )
    .unwrap();
    db
}

#[test]
fn fastq_fai_seek_returns_correct_record() {
    // read3 is the 3rd record — a seek should land on it directly, not stream from the top
    assert_eq!(
        scalar_str(
            &fq_fai_db(),
            "SELECT sequence FROM fq_fai WHERE id = 'read3'"
        ),
        "ACGTTTTT"
    );
}

#[test]
fn fastq_fai_seek_does_not_bleed_into_adjacent_record() {
    let db = fq_fai_db();
    // Only one row should come back — confirms the reader stops after the seeked record
    assert_eq!(
        scalar_i64(&db, "SELECT COUNT(*) FROM fq_fai WHERE id = 'read2'"),
        1
    );
    assert_eq!(
        scalar_str(&db, "SELECT sequence FROM fq_fai WHERE id = 'read2'"),
        "GGGG"
    );
}

#[test]
fn fastq_fai_seek_first_record() {
    // read1 is at byte 0 — the backward scan in find_record_offset must handle hitting start-of-file
    let db = fq_fai_db();
    assert_eq!(
        scalar_str(&db, "SELECT sequence FROM fq_fai WHERE id = 'read1'"),
        "ACGT"
    );
    assert_eq!(
        scalar_i64(&db, "SELECT COUNT(*) FROM fq_fai WHERE id = 'read1'"),
        1
    );
}

#[test]
fn fastq_fai_seek_last_record_preserves_lowercase() {
    // read4 is the last record and has lowercase bases — tests boundary + case preservation
    assert_eq!(
        scalar_str(
            &fq_fai_db(),
            "SELECT sequence FROM fq_fai WHERE id = 'read4'"
        ),
        "acgtgggg"
    );
}

#[test]
fn fastq_fai_seek_quality_correct() {
    // Quality string must match the seeked record exactly, not a neighbouring record's quality
    let db = fq_fai_db();
    assert_eq!(
        scalar_str(&db, "SELECT quality FROM fq_fai WHERE id = 'read3'"),
        "????????"
    );
    assert_eq!(
        scalar_str(&db, "SELECT quality FROM fq_fai WHERE id = 'read1'"),
        "IIII"
    );
}

#[test]
fn fastq_fai_seek_unknown_id_returns_empty() {
    // A miss in the FAI must fall back to a full scan and return no rows, not an error
    assert_eq!(
        scalar_i64(
            &fq_fai_db(),
            "SELECT COUNT(*) FROM fq_fai WHERE id = 'nonexistent'"
        ),
        0
    );
}

#[test]
fn fastq_fai_sequential_seeks_reset_correctly() {
    // Two id=Eq queries on the same table reuse the cursor; exit_early must be cleared between
    // them or the second query immediately returns EOF.
    let db = fq_fai_db();
    assert_eq!(
        scalar_str(&db, "SELECT sequence FROM fq_fai WHERE id = 'read1'"),
        "ACGT"
    );
    assert_eq!(
        scalar_str(&db, "SELECT sequence FROM fq_fai WHERE id = 'read4'"),
        "acgtgggg"
    );
    assert_eq!(
        scalar_str(&db, "SELECT sequence FROM fq_fai WHERE id = 'read2'"),
        "GGGG"
    );
}

// --- has_stop_codon scalar ---

#[test]
fn has_stop_codon_dna_taa() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('TAA')"), 1);
}

#[test]
fn has_stop_codon_dna_tag() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('TAG')"), 1);
}

#[test]
fn has_stop_codon_dna_tga() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('TGA')"), 1);
}

#[test]
fn has_stop_codon_rna_uaa() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('UAA')"), 1);
}

#[test]
fn has_stop_codon_rna_uag() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('UAG')"), 1);
}

#[test]
fn has_stop_codon_rna_uga() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('UGA')"), 1);
}

#[test]
fn has_stop_codon_lowercase() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('taa')"), 1);
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('tag')"), 1);
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('tga')"), 1);
}

#[test]
fn has_stop_codon_embedded() {
    // stop codon not at position 0
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('ACGTAA')"), 1);
}

#[test]
fn has_stop_codon_no_stop() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('ACGT')"), 0);
}

#[test]
fn has_stop_codon_empty() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('')"), 0);
}

#[test]
fn has_stop_codon_too_short() {
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('TA')"), 0);
}

#[test]
fn has_stop_codon_near_miss() {
    // TGG is Trp, not a stop codon
    assert_eq!(scalar_i64(&db(), "SELECT has_stop_codon('ACGTGG')"), 0);
}

#[test]
fn has_stop_codon_on_fasta_table() {
    // none of the fixture sequences contain a stop codon
    assert_eq!(
        scalar_i64(
            &fasta_db(),
            "SELECT COUNT(*) FROM fa WHERE has_stop_codon(sequence) = 1"
        ),
        0
    );
}

// --- bgzf FASTA: streaming + future seek (Pulchroboletus rubricitrinus ITS, MF193883-MF193886) ---
// Sequences from: Farid, Franck & Garey (2017), Czech Mycology 69(2):143-162.
// Fixture is bgzf-compressed with a .fai alongside.
// Currently exercises the streaming path; seek path covered once bgzf seek is implemented.

fn bgzf_db() -> Database {
    assert!(
        std::path::Path::new(TEST_BGZF_FA).exists(),
        "bgzf fixture missing: {TEST_BGZF_FA}"
    );
    assert!(
        std::path::Path::new(TEST_BGZF_FAI).exists(),
        "bgzf fai fixture missing: {TEST_BGZF_FAI}"
    );
    let db = db();
    db.execute(
        &format!("CREATE VIRTUAL TABLE bgzf USING fasta('{TEST_BGZF_FA}')"),
        (),
    )
    .unwrap();
    db
}

#[test]
fn bgzf_row_count() {
    assert_eq!(scalar_i64(&bgzf_db(), "SELECT COUNT(*) FROM bgzf"), 4);
}

#[test]
fn bgzf_first_record_length() {
    assert_eq!(
        scalar_i64(
            &bgzf_db(),
            "SELECT length FROM bgzf WHERE id = 'MF193883.1'"
        ),
        715
    );
}

#[test]
fn bgzf_last_record_length() {
    assert_eq!(
        scalar_i64(
            &bgzf_db(),
            "SELECT length FROM bgzf WHERE id = 'MF193886.1'"
        ),
        732
    );
}

#[test]
fn bgzf_record_with_ambiguous_base() {
    // MF193884.1 contains an N at position 4
    assert_eq!(
        scalar_i64(
            &bgzf_db(),
            "SELECT length FROM bgzf WHERE id = 'MF193884.1'"
        ),
        721
    );
}

#[test]
fn bgzf_conserved_its_region_matches_all() {
    // All four ITS sequences share this conserved 5.8S/ITS2 boundary region
    assert_eq!(
        scalar_i64(
            &bgzf_db(),
            "SELECT COUNT(*) FROM bgzf WHERE sequence LIKE '%GCATCGATGAAGAACGCAGCGAATTGCGAT%'"
        ),
        4
    );
}

#[test]
fn bgzf_unknown_id_returns_empty() {
    assert_eq!(
        scalar_i64(
            &bgzf_db(),
            "SELECT COUNT(*) FROM bgzf WHERE id = 'nonexistent'"
        ),
        0
    );
}

#[test]
fn bgzf_sequence_no_match_returns_empty() {
    assert_eq!(
        scalar_i64(
            &bgzf_db(),
            "SELECT COUNT(*) FROM bgzf WHERE sequence LIKE '%AAAAAAAAAAAAAAAAAAAAAA%'"
        ),
        0
    );
}

#[test]
fn bgzf_fai_seek_middle_record_correct_sequence() {
    // MF193884.1 is the second record — seeks past the first block, verifies correct landing
    let seq = scalar_str(
        &bgzf_db(),
        "SELECT sequence FROM bgzf WHERE id = 'MF193884.1'",
    );
    assert!(
        seq.starts_with("TTCNGTAGGG"),
        "unexpected sequence start: {seq}"
    );
    assert_eq!(seq.len(), 721);
}

#[test]
fn bgzf_fai_seek_does_not_bleed_into_adjacent_record() {
    // Seeking to a middle record must return exactly one row with the correct id
    let count = scalar_i64(
        &bgzf_db(),
        "SELECT COUNT(*) FROM bgzf WHERE id = 'MF193885.1'",
    );
    assert_eq!(count, 1);
    let seq = scalar_str(
        &bgzf_db(),
        "SELECT sequence FROM bgzf WHERE id = 'MF193885.1'",
    );
    assert_eq!(seq.len(), 726);
    assert!(
        seq.starts_with("CTTCCGTAGGT"),
        "unexpected sequence start: {seq}"
    );
}

// --- bgzf FASTQ: seek + edge cases ---
// Fixture: sample.fastq.gz (4 records, all length 12)
//   read1  ACGTACGTACGT  quality IIIIIIIIIIII (mean 40, min 40)
//   read2  GGTTGGTTGGTT  quality !!!!!!!!!!!!! (mean 0, min 0)
//   read3  AAAAACCCCCGG  quality ????????IIII (mean ~33, min 30)
//   read4  acgtacgtacgt  quality IIIIIIIIIIII (lowercase sequence)

fn bgzf_fastq_db() -> Database {
    assert!(
        std::path::Path::new(TEST_BGZF_FASTQ).exists(),
        "bgzf fastq fixture missing: {TEST_BGZF_FASTQ}"
    );
    assert!(
        std::path::Path::new(TEST_BGZF_FASTQ_FAI).exists(),
        "bgzf fastq fai fixture missing: {TEST_BGZF_FASTQ_FAI}"
    );
    let db = db();
    db.execute(
        &format!("CREATE VIRTUAL TABLE bgzf_fq USING fastq('{TEST_BGZF_FASTQ}')"),
        (),
    )
    .unwrap();
    db
}

#[test]
fn bgzf_fastq_row_count() {
    assert_eq!(
        scalar_i64(&bgzf_fastq_db(), "SELECT COUNT(*) FROM bgzf_fq"),
        4
    );
}

#[test]
fn bgzf_fastq_seek_first_record() {
    assert_eq!(
        scalar_str(
            &bgzf_fastq_db(),
            "SELECT sequence FROM bgzf_fq WHERE id = 'read1'"
        ),
        "ACGTACGTACGT"
    );
}

#[test]
fn bgzf_fastq_seek_middle_record() {
    assert_eq!(
        scalar_str(
            &bgzf_fastq_db(),
            "SELECT sequence FROM bgzf_fq WHERE id = 'read2'"
        ),
        "GGTTGGTTGGTT"
    );
}

#[test]
fn bgzf_fastq_seek_last_record() {
    assert_eq!(
        scalar_str(
            &bgzf_fastq_db(),
            "SELECT sequence FROM bgzf_fq WHERE id = 'read4'"
        ),
        "acgtacgtacgt"
    );
}

#[test]
fn bgzf_fastq_seek_does_not_bleed_into_adjacent() {
    let count = scalar_i64(
        &bgzf_fastq_db(),
        "SELECT COUNT(*) FROM bgzf_fq WHERE id = 'read3'",
    );
    assert_eq!(count, 1);
    assert_eq!(
        scalar_str(
            &bgzf_fastq_db(),
            "SELECT sequence FROM bgzf_fq WHERE id = 'read3'"
        ),
        "AAAAACCCCCGG"
    );
}

#[test]
fn bgzf_fastq_quality_preserved_after_seek() {
    assert_eq!(
        scalar_str(
            &bgzf_fastq_db(),
            "SELECT quality FROM bgzf_fq WHERE id = 'read1'"
        ),
        "IIIIIIIIIIII"
    );
}

#[test]
fn bgzf_fastq_seek_unknown_id_returns_empty() {
    assert_eq!(
        scalar_i64(
            &bgzf_fastq_db(),
            "SELECT COUNT(*) FROM bgzf_fq WHERE id = 'nonexistent'"
        ),
        0
    );
}

#[test]
fn bgzf_fastq_quality_filter_streaming() {
    // mean_quality > 30 matches read1 (40), read3 (~33), read4 (40); excludes read2 (0)
    assert_eq!(
        scalar_i64(
            &bgzf_fastq_db(),
            "SELECT COUNT(*) FROM bgzf_fq WHERE mean_quality > 30"
        ),
        3
    );
}

#[test]
fn bgzf_fastq_zero_quality_record() {
    // read2 has all '!' quality (Phred 0)
    assert_eq!(
        scalar_i64(
            &bgzf_fastq_db(),
            "SELECT COUNT(*) FROM bgzf_fq WHERE min_quality = 0"
        ),
        1
    );
    assert_eq!(
        scalar_str(
            &bgzf_fastq_db(),
            "SELECT id FROM bgzf_fq WHERE min_quality = 0"
        ),
        "read2"
    );
}

#[test]
fn bgzf_fastq_lowercase_sequence_preserved_after_seek() {
    // read4 has lowercase sequence — seek must not normalise it
    let seq = scalar_str(
        &bgzf_fastq_db(),
        "SELECT sequence FROM bgzf_fq WHERE id = 'read4'",
    );
    assert_eq!(seq, "acgtacgtacgt");
}

// --- n50 aggregate: empty result set ---

#[test]
fn n50_empty_result_set() {
    assert_eq!(
        scalar_i64(&db(), "SELECT n50(x) FROM (SELECT 1 AS x) WHERE x < 0"),
        0
    );
}
