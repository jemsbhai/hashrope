//! Registry-only release smoke tests. No repository or path dependencies.

#[cfg(test)]
mod tests {
    use hashrope::{Arena, PolynomialHash, MERSENNE_61};
    use hashrope_bio::gene_diff::{
        diff_genes, diff_genes_single_arena, screen_genes_single_arena, Exon, GeneRegion,
    };
    use hashrope_bio::resistance::{
        check_panel_single_arena, check_resistance_panel, ResistanceSite,
    };
    use std::panic::{catch_unwind, AssertUnwindSafe};

    fn gene() -> GeneRegion {
        GeneRegion {
            name: "test".into(),
            transcript: "test".into(),
            chrom: "chr1".into(),
            strand: '+',
            tx_start: 0,
            tx_end: 3,
            exons: vec![Exon {
                index: 1,
                start: 0,
                end: 3,
            }],
        }
    }

    fn site() -> ResistanceSite {
        ResistanceSite {
            gene: "test",
            position: 1,
            wildtype_aa: "X",
            mutant_aa: "Y",
            drug_class: "test",
            annotation: "byte regression",
        }
    }

    #[test]
    fn default_byte_and_token_collisions_are_exact() {
        for (a, b, fingerprint) in [
            (vec![0, 131], vec![1, 0], 263),
            (
                33536u32.to_le_bytes().to_vec(),
                1u32.to_le_bytes().to_vec(),
                4_513_475,
            ),
        ] {
            let mut arena = Arena::new();
            let left = arena.from_bytes(&a);
            let right = arena.from_bytes(&b);
            assert_eq!(arena.hash(left), fingerprint);
            assert_eq!(arena.hash(right), fingerprint);
            assert!(!hashrope::rope_equal(&arena, left, right));
            assert_eq!(hashrope::rope_lcp(&arena, left, right), 0);
            assert!(!arena.equal_range(left, 0, right, 0, a.len() as u64));
        }
    }

    #[test]
    fn cross_profile_ranges_repeats_and_bounds() {
        let mut a = Arena::new_lazy();
        let pattern = a.from_bytes(b"ab");
        let repeated = a.repeat(pattern, 1 << 40);
        let count = a.node_count();
        let mut b = Arena::with_hash(MERSENNE_61, 257);
        let right = b.from_bytes(b"zababx");
        assert!(a.equal_range_across(repeated, 0, &b, right, 1, 4));
        assert_eq!(a.lcp_range_across(repeated, 0, &b, right, 1, 5), 4);
        assert!(a.equal_range(repeated, 1, repeated, 3, (1 << 40) - 3));
        assert_eq!(a.node_count(), count);
        assert!(catch_unwind(|| a.equal_range(repeated, u64::MAX, repeated, 0, 2)).is_err());
        assert!(catch_unwind(|| a.equal_range(None, 1, None, 1, 0)).is_err());
    }

    #[test]
    fn changed_profile_is_rejected_without_blocking_exact_bytes() {
        let mut arena = Arena::new();
        let a = arena.from_bytes(b"abc");
        let b = arena.from_bytes(b"abc");
        let fingerprint = arena.hash(a);
        let count = arena.node_count();
        *arena.hasher_mut() = PolynomialHash::new(MERSENNE_61, 257);
        assert!(arena.equal(a, b));
        assert_eq!(arena.lcp(a, b), 3);
        assert!(catch_unwind(|| arena.hash(a)).is_err());
        assert!(catch_unwind(AssertUnwindSafe(|| arena.concat(a, b))).is_err());
        assert_eq!(arena.node_count(), count);
        *arena.hasher_mut() = PolynomialHash::new(MERSENNE_61, 131);
        assert_eq!(arena.hash(a), fingerprint);
    }

    #[test]
    fn biology_detects_colliding_gene_exon_and_codon() {
        let mut arena = Arena::new();
        let a = arena.from_bytes(&[0, 131, 42]);
        let b = arena.from_bytes(&[1, 0, 42]);
        assert_eq!(arena.hash(a), arena.hash(b));
        let report = diff_genes_single_arena(&mut arena, a, b, &[gene()], None, true);
        assert!(report.gene_diffs[0].body_changed());
        assert!(report.gene_diffs[0].exon_diffs[0].changed());
        assert_eq!(
            screen_genes_single_arena(&mut arena, a, b, &[gene()], None),
            vec!["test"]
        );
        let panel = check_panel_single_arena(&mut arena, a, b, &[site()]);
        assert_eq!(panel[0].ref_hash, panel[0].sample_hash);
        assert!(panel[0].is_mutant);
    }

    #[test]
    fn biology_ignores_profile_differences_for_equal_content() {
        let mut a = Arena::new_lazy();
        let mut b = Arena::with_hash(MERSENNE_61, 257);
        let left = a.from_bytes(b"ACG");
        let right = b.from_bytes(b"ACG");
        let report = diff_genes(&mut a, left, &mut b, right, &[gene()], None);
        assert_ne!(
            report.gene_diffs[0].body_ref_hash,
            report.gene_diffs[0].body_sample_hash
        );
        assert!(!report.gene_diffs[0].body_changed());
        assert!(!report.gene_diffs[0].exon_diffs[0].changed());
        assert!(!check_resistance_panel(&mut a, left, &mut b, right, &[site()])[0].is_mutant);
    }

    #[cfg(feature = "bigint")]
    #[test]
    fn bigint_huge_repeats_ranges_and_profile_guards() {
        use hashrope::bigint::{
            Arena as BigArena, BigUint, PolynomialHash as BigHash, MERSENNE_127,
        };
        let n = |value: u32| BigUint::from(value);
        let copies = n(1) << 1300usize;
        let mut a = BigArena::with_hash(n(127), n(31));
        let mut b = BigArena::with_hash(MERSENNE_127.into(), n(257));
        let ab = a.from_bytes(b"ab");
        let abab = b.from_bytes_chunked(b"abab", 1);
        let x = a.repeat(ab, &copies);
        let y = b.repeat(abab, &(&copies >> 1usize));
        let count = (a.node_count(), b.node_count());
        assert!(a.equal_across(x, &b, y));
        assert_eq!(a.lcp_across(x, &b, y), &copies * 2u8);
        let start = (n(1) << 1100usize) + 1u8;
        assert!(a.equal_range_across(x, &start, &b, y, &(&start + 2u8), &(&copies - 3u8)));
        assert_eq!((a.node_count(), b.node_count()), count);
        assert!(catch_unwind(|| a.equal_range(x, &(&copies * 3u8), x, &n(0), &n(0))).is_err());
        *a.hasher_mut() = BigHash::new(n(127), n(37));
        assert!(a.equal_across(x, &b, y));
        assert!(catch_unwind(|| a.hash(x)).is_err());
        *a.hasher_mut() = BigHash::new(n(127), n(31));
        assert_eq!(a.len(x), &copies * 2u8);
    }

    #[cfg(feature = "bigint")]
    #[test]
    fn bigint_default_fingerprint_collision_is_not_equality() {
        let mut arena = hashrope::bigint::Arena::new();
        let a = arena.from_bytes(&[0, 131]);
        let b = arena.from_bytes(&[1, 0]);
        assert_eq!(arena.hash(a), 263u32.into());
        assert_eq!(arena.hash(a), arena.hash(b));
        assert!(!arena.equal(a, b));
        assert_eq!(arena.lcp(a, b), 0u32.into());
    }
}
