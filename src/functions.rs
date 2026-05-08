pub fn compute_gc(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
        .count();
    gc as f64 / seq.len() as f64
}

pub fn dna_to_rna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'T' => b'U',
            b't' => b'u',
            _ => b,
        })
        .collect()
}

pub fn rna_to_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'U' => b'T',
            b'u' => b't',
            _ => b,
        })
        .collect()
}

pub fn n_count(seq: &[u8]) -> i64 {
    seq.iter().filter(|&b| matches!(b, b'n' | b'N')).count() as i64
}

pub fn base_count(seq: &[u8], base: u8) -> sqlite3_ext::Result<i64> {
    if !matches!(
        base.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T' | b'U' | b'N'
    ) {
        return Err("Not a valid base".into());
    }
    let base = base.to_ascii_uppercase();
    Ok(seq
        .iter()
        .filter(|b| b.to_ascii_uppercase() == base)
        .count() as i64)
}

pub fn is_valid_dna(seq: &[u8]) -> bool {
    seq.iter()
        .all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}
pub fn is_valid_rna(seq: &[u8]) -> bool {
    seq.iter()
        .all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'U' | b'N'))
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'a' => b't',
            b'G' => b'C',
            b'g' => b'c',
            b'T' => b'A',
            b't' => b'a',
            b'C' => b'G',
            b'c' => b'g',
            b'U' => b'A',
            b'u' => b'a',
            _ => b,
        })
        .collect()
}

pub fn mean_quality(qual: &[u8]) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum: f64 = qual.iter().map(|&b| (b - 33) as f64).sum();
    sum / qual.len() as f64
}

pub fn min_quality(qual: &[u8]) -> i64 {
    qual.iter().map(|&b| (b - 33) as i64).min().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    // compute_gc tests
    #[test]
    fn gc_empty() {
        assert_eq!(compute_gc(b""), 0.0);
    }

    #[test]
    fn gc_all_gc() {
        assert_eq!(compute_gc(b"GGCC"), 1.0);
        assert_eq!(compute_gc(b"ggcc"), 1.0);
    }

    #[test]
    fn gc_no_gc() {
        assert_eq!(compute_gc(b"AATT"), 0.0);
        assert_eq!(compute_gc(b"aatt"), 0.0);
    }

    #[test]
    fn gc() {
        assert_eq!(compute_gc(b"ACGT"), 0.5);
    }

    #[test]
    fn gc_with_ambiguous_bases() {
        assert_eq!(compute_gc(b"ATGCN"), 0.4);
    }

    // n_count
    #[test]
    fn ncount_empty() {
        assert_eq!(n_count(b""), 0);
    }

    #[test]
    fn ncount_all_n() {
        assert_eq!(n_count(b"NNNN"), 4);
    }

    #[test]
    fn ncount_no_n() {
        assert_eq!(n_count(b"ACGT"), 0);
    }

    #[test]
    fn ncount_3n() {
        assert_eq!(n_count(b"ACNGNTN"), 3);
    }

    #[test]
    fn ncount_lower_n() {
        assert_eq!(n_count(b"ACnGnTn"), 3);
    }

    // base_count
    #[test]
    fn base_count_basic() {
        assert_eq!(base_count(b"ACGT", b'A'), Ok(1));
        assert_eq!(base_count(b"AAAA", b'A'), Ok(4));
    }

    #[test]
    fn base_count_case_insensitive() {
        assert_eq!(base_count(b"AaCcGgTt", b'a'), Ok(2));
        assert_eq!(base_count(b"AaCcGgTt", b'A'), Ok(2));
    }

    #[test]
    fn base_count_invalid_base() {
        assert!(base_count(b"ACGT", b'Z').is_err());
    }

    #[test]
    fn base_count_empty() {
        assert_eq!(base_count(b"", b'A'), Ok(0));
    }

    //to_rna
    #[test]
    fn rna_converts_t_to_u() {
        assert_eq!(dna_to_rna(b"ACGT"), b"ACGU");
    }

    #[test]
    fn rna_lowercase_t() {
        assert_eq!(dna_to_rna(b"acgt"), b"acgu");
    }

    #[test]
    fn rna_no_t() {
        assert_eq!(dna_to_rna(b"ACGA"), b"ACGA");
    }

    #[test]
    fn rna_passthrough_non_dna() {
        assert_eq!(dna_to_rna(b"AFCGA"), b"AFCGA");
    }

    //Reverse Complement
    #[test]
    fn reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn reverse_complement_single_bases() {
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"T"), b"A");
        assert_eq!(reverse_complement(b"G"), b"C");
        assert_eq!(reverse_complement(b"C"), b"G");
        assert_eq!(reverse_complement(b"U"), b"A");
    }

    #[test]
    fn reverse_complement_ambiguous_passthrough() {
        assert_eq!(reverse_complement(b"N"), b"N");
        assert_eq!(reverse_complement(b"n"), b"n");
    }

    #[test]
    fn reverse_complement_idempotent() {
        let seq = b"ACGTACGT";
        assert_eq!(reverse_complement(&reverse_complement(seq)), seq);
    }

    #[test]
    fn reverse_complement_pure_dna() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
    }

    #[test]
    fn reverse() {
        assert_eq!(reverse_complement(b"AGCTUagctuNn"), b"nNaagctAAGCT")
    }

    // is_valid_dna
    #[test]
    fn valid_dna_basic() {
        assert!(is_valid_dna(b"ACGT"));
    }

    #[test]
    fn valid_dna_lowercase() {
        assert!(is_valid_dna(b"acgt"));
    }

    #[test]
    fn valid_dna_mixed_case() {
        assert!(is_valid_dna(b"AcGt"));
    }

    #[test]
    fn valid_dna_with_n() {
        assert!(is_valid_dna(b"ACGTN"));
        assert!(is_valid_dna(b"acgtn"));
    }

    #[test]
    fn valid_dna_empty() {
        assert!(is_valid_dna(b""));
    }

    #[test]
    fn valid_dna_rejects_u() {
        assert!(!is_valid_dna(b"ACGTU"));
    }

    #[test]
    fn valid_dna_rejects_invalid() {
        assert!(!is_valid_dna(b"ACGTZ"));
        assert!(!is_valid_dna(b"ACGT1"));
        assert!(!is_valid_dna(b"ACGT!"));
    }

    #[test]
    fn valid_dna_rejects_on_first_invalid() {
        // invalid base at start — should short circuit
        assert!(!is_valid_dna(b"ZACGT"));
    }

    #[test]
    fn valid_dna_single_base() {
        assert!(is_valid_dna(b"A"));
        assert!(is_valid_dna(b"C"));
        assert!(is_valid_dna(b"G"));
        assert!(is_valid_dna(b"T"));
        assert!(is_valid_dna(b"N"));
        assert!(!is_valid_dna(b"U"));
        assert!(!is_valid_dna(b"Z"));
    }

    // is_valid_rna
    #[test]
    fn valid_rna_basic() {
        assert!(is_valid_rna(b"ACGU"));
    }

    #[test]
    fn valid_rna_lowercase() {
        assert!(is_valid_rna(b"acgu"));
    }

    #[test]
    fn valid_rna_mixed_case() {
        assert!(is_valid_rna(b"AcGu"));
    }

    #[test]
    fn valid_rna_with_n() {
        assert!(is_valid_rna(b"ACGUN"));
        assert!(is_valid_rna(b"acgun"));
    }

    #[test]
    fn valid_rna_empty() {
        assert!(is_valid_rna(b""));
    }

    #[test]
    fn valid_rna_rejects_t() {
        assert!(!is_valid_rna(b"ACGT"));
    }

    #[test]
    fn valid_rna_rejects_invalid() {
        assert!(!is_valid_rna(b"ACGUZ"));
        assert!(!is_valid_rna(b"ACGU1"));
        assert!(!is_valid_rna(b"ACGU!"));
    }

    #[test]
    fn valid_rna_single_base() {
        assert!(is_valid_rna(b"A"));
        assert!(is_valid_rna(b"C"));
        assert!(is_valid_rna(b"G"));
        assert!(is_valid_rna(b"U"));
        assert!(is_valid_rna(b"N"));
        assert!(!is_valid_rna(b"T"));
        assert!(!is_valid_rna(b"Z"));
    }

    #[test]
    fn mean_quality_basic() {
        // All 'I' = ASCII 73, Phred 40
        assert_eq!(mean_quality(b"IIII"), 40.0);
    }

    #[test]
    fn mean_quality_min_score() {
        // All '!' = ASCII 33, Phred 0
        assert_eq!(mean_quality(b"!!!!"), 0.0);
    }

    #[test]
    fn mean_quality_mixed() {
        // '!' = Q0, 'I' = Q40, mean = 20.0
        assert_eq!(mean_quality(b"!I"), 20.0);
    }

    #[test]
    fn mean_quality_single_base() {
        // '5' = ASCII 53, Phred 20
        assert_eq!(mean_quality(b"5"), 20.0);
    }

    #[test]
    fn mean_quality_q30() {
        // '?' = ASCII 63, Phred 30
        assert_eq!(mean_quality(b"????"), 30.0);
    }

    #[test]
    fn min_quality_basic() {
        // 'I' = Q40, '!' = Q0
        assert_eq!(min_quality(b"III!III"), 0);
    }

    #[test]
    fn min_quality_all_same() {
        assert_eq!(min_quality(b"????"), 30);
    }

    #[test]
    fn min_quality_single_base() {
        // '5' = ASCII 53, Phred 20
        assert_eq!(min_quality(b"5"), 20);
    }

    #[test]
    fn min_quality_returns_lowest() {
        // '+' = ASCII 43, Phred 10 — lowest in mixed string
        assert_eq!(min_quality(b"IIII+IIII"), 10);
    }

    #[test]
    fn mean_quality_empty() {
        assert_eq!(mean_quality(b""), 0.0);
    }

    #[test]
    fn min_quality_empty() {
        assert_eq!(min_quality(b""), 0);
    }
}
