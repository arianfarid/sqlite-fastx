// Scalar

use std::fmt::Display;

use sqlite3_ext::{
    FromValue,
    function::{AggregateFunction, FromUserData, ToContextResult},
};

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

pub struct BaseComposition {
    pub a: f64,
    pub c: f64,
    pub g: f64,
    pub t: f64,
    pub u: f64,
}
impl Display for BaseComposition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{")?;
        let mut sep = "";
        for (base, val) in [
            ("A", self.a),
            ("C", self.c),
            ("G", self.g),
            ("T", self.t),
            ("U", self.u),
        ] {
            write!(f, "{}\"{}\": {:.4}", sep, base, val)?;
            sep = ", ";
        }
        write!(f, "}}")
    }
}
pub fn base_composition(seq: &[u8]) -> BaseComposition {
    let mut bcomp = BaseComposition {
        a: 0.,
        c: 0.,
        g: 0.,
        t: 0.,
        u: 0.,
    };
    if seq.is_empty() {
        return bcomp;
    }
    let total = seq.len() as f64;
    seq.iter().for_each(|&b| match &b.to_ascii_uppercase() {
        b'A' => bcomp.a += 1.,
        b'G' => bcomp.g += 1.,
        b'T' => bcomp.t += 1.,
        b'C' => bcomp.c += 1.,
        b'U' => bcomp.u += 1.,
        _ => {}
    });

    bcomp.a /= total;
    bcomp.g /= total;
    bcomp.c /= total;
    bcomp.t /= total;
    bcomp.u /= total;

    bcomp
}

// Aggregates

pub struct N50Accumulator {
    lengths: Vec<i64>,
}
impl FromUserData<()> for N50Accumulator {
    fn from_user_data(_: &()) -> Self {
        N50Accumulator { lengths: vec![] }
    }
}
impl AggregateFunction<()> for N50Accumulator {
    fn step(
        &mut self,
        _context: &sqlite3_ext::function::Context,
        args: &mut [&mut sqlite3_ext::ValueRef],
    ) -> sqlite3_ext::Result<()> {
        let value = args[0].get_i64();
        self.lengths.push(value);
        Ok(())
    }

    fn value(&self, context: &sqlite3_ext::function::Context) -> sqlite3_ext::Result<()> {
        if self.lengths.is_empty() {
            return context.set_result(0i64);
        }
        let mut sorted = self.lengths.clone();
        sorted.sort_unstable_by(|a, b| b.cmp(a));
        let total: i64 = self.lengths.clone().iter().sum();
        let half = total / 2;
        let mut running = 0i64;
        for &len in &sorted {
            running += len;
            if running >= half {
                return context.set_result(len);
            }
        }
        context.set_result(0i64)
    }

    fn inverse(
        &mut self,
        _context: &sqlite3_ext::function::Context,
        args: &mut [&mut sqlite3_ext::ValueRef],
    ) -> sqlite3_ext::Result<()> {
        let value = args[0].get_i64();
        let pos = self.lengths.iter().position(|&x| x == value);
        match pos {
            Some(position) => {
                self.lengths.remove(position);
            }
            None => {}
        };
        Ok(())
    }

    fn default_value(
        user_data: &(),
        context: &sqlite3_ext::function::Context,
    ) -> sqlite3_ext::Result<()>
    where
        Self: Sized,
    {
        Self::from_user_data(user_data).value(context)
    }
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

    // base_composition
    #[test]
    fn base_composition_empty() {
        let bc = base_composition(b"");
        assert_eq!(bc.a, 0.0);
        assert_eq!(bc.c, 0.0);
        assert_eq!(bc.g, 0.0);
        assert_eq!(bc.t, 0.0);
        assert_eq!(bc.u, 0.0);
    }

    #[test]
    fn base_composition_equal_dna() {
        let bc = base_composition(b"ACGT");
        assert_eq!(bc.a, 0.25);
        assert_eq!(bc.c, 0.25);
        assert_eq!(bc.g, 0.25);
        assert_eq!(bc.t, 0.25);
        assert_eq!(bc.u, 0.0);
    }

    #[test]
    fn base_composition_equal_rna() {
        let bc = base_composition(b"ACGU");
        assert_eq!(bc.a, 0.25);
        assert_eq!(bc.c, 0.25);
        assert_eq!(bc.g, 0.25);
        assert_eq!(bc.t, 0.0);
        assert_eq!(bc.u, 0.25);
    }

    #[test]
    fn base_composition_all_one_base() {
        let bc = base_composition(b"AAAA");
        assert_eq!(bc.a, 1.0);
        assert_eq!(bc.c, 0.0);
        assert_eq!(bc.g, 0.0);
        assert_eq!(bc.t, 0.0);
        assert_eq!(bc.u, 0.0);
    }

    #[test]
    fn base_composition_lowercase() {
        let bc = base_composition(b"acgt");
        assert_eq!(bc.a, 0.25);
        assert_eq!(bc.c, 0.25);
        assert_eq!(bc.g, 0.25);
        assert_eq!(bc.t, 0.25);
        assert_eq!(bc.u, 0.0);
    }

    #[test]
    fn base_composition_mixed_case() {
        let upper = base_composition(b"AACCGGTT");
        let lower = base_composition(b"aaccggtt");
        assert_eq!(upper.a, lower.a);
        assert_eq!(upper.c, lower.c);
        assert_eq!(upper.g, lower.g);
        assert_eq!(upper.t, lower.t);
    }

    #[test]
    fn base_composition_with_n() {
        // N is not counted in any base but is included in total length,
        // so fractions are diluted
        let bc = base_composition(b"ACGTN");
        assert_eq!(bc.a, 0.2);
        assert_eq!(bc.c, 0.2);
        assert_eq!(bc.g, 0.2);
        assert_eq!(bc.t, 0.2);
        assert_eq!(bc.u, 0.0);
    }

    #[test]
    fn base_composition_sums_to_one_dna() {
        let bc = base_composition(b"AACCGGTT");
        let sum = bc.a + bc.c + bc.g + bc.t + bc.u;
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn base_composition_sums_to_one_rna() {
        let bc = base_composition(b"AACCGGUU");
        let sum = bc.a + bc.c + bc.g + bc.t + bc.u;
        assert!((sum - 1.0).abs() < 1e-10);
    }
}
