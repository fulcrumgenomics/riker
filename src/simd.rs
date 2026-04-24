//! Reusable SIMD kernels for riker's hot byte-level loops, built on the [`wide`] crate.
//!
//! Each kernel processes a fixed-width chunk of its input in SIMD and handles the
//! sub-chunk tail with a plain scalar loop. The kernels are intentionally general
//! — "count bases ≥ Q over a slice", "GC count of a slice" — rather than tailored
//! to any one command, so a new call site can adopt them without growing the module.
//!
//! ## Conventions
//!
//! * **Raw Phred quality scores.** noodles decodes BAM/SAM/CRAM qualities to
//!   raw Phred bytes, so quality-threshold kernels compare directly
//!   (`q >= threshold`) — there is no `+33` offset to strip.
//! * **Case-folding via `| 0x20`.** Bit 5 distinguishes ASCII upper- and lower-case
//!   letters, so `byte | 0x20` collapses `A..Z` onto `a..z`. Non-letter bytes keep
//!   their identity under this fold (the `| 0x20` may flip a non-letter bit, but
//!   the subsequent equality test against a lowercase letter still fails).
//! * **Lane width.** All SIMD kernels here process per-read byte slices
//!   (typical lengths 100–300 bytes) and use `u8x16` to keep the scalar tail
//!   short. Apple Silicon NEON is 128-bit so `u8x16` is native; on AVX2 x86 it
//!   maps to XMM registers. Widening to `u8x32` did not pay for itself at
//!   these slice sizes during benchmarking.
//! * **Empty slices** return the neutral value (`0` for counts).

// The `chunks_exact(N)` → `chunk.try_into::<[u8; N]>().unwrap()` pattern is
// infallible by construction — `chunks_exact` guarantees every yielded chunk
// has length `N`. Silence the pedantic `missing_panics_doc` lint at the
// module level (including tests) rather than adding empty `# Panics`
// sections to every kernel.
#![allow(clippy::missing_panics_doc)]

use wide::u8x16;

/// Count quality bytes in `qual` with value `>= threshold`.
///
/// Qualities are raw Phred values, not Phred+33.
#[must_use]
pub fn count_bases_ge_q(qual: &[u8], threshold: u8) -> u64 {
    let cutoff = u8x16::splat(threshold);
    let mut count = 0u64;
    let chunks = qual.chunks_exact(16);
    let tail = chunks.remainder();
    for chunk in chunks {
        let v = u8x16::new(chunk.try_into().unwrap());
        count += u64::from(v.simd_ge(cutoff).to_bitmask().count_ones());
    }
    for &q in tail {
        if q >= threshold {
            count += 1;
        }
    }
    count
}

/// Count quality bytes in `qual` with value `< threshold`.
///
/// Dual of [`count_bases_ge_q`] — useful for low-quality filters where the
/// excluded count is the natural output (e.g. `bases_excl_baseq`). Expressed in
/// terms of [`count_bases_ge_q`] because the two partitions of the slice must
/// always sum to `qual.len()`; no separate SIMD kernel is needed.
#[must_use]
// Inline so the sole `sub` folds into the caller's expression and this
// wrapper disappears from the call graph.
#[inline]
pub fn count_bases_lt_q(qual: &[u8], threshold: u8) -> u64 {
    qual.len() as u64 - count_bases_ge_q(qual, threshold)
}

/// Count G and C bases in `seq`, matching both upper- and lower-case.
///
/// Non-ACGT bytes do not match. Ambiguity codes other than `N` (e.g. `S`, which
/// encodes "G or C") are treated as non-GC, matching the behavior of the scalar
/// `match b { b'G' | b'C' | b'g' | b'c' => … }` form it replaces.
#[must_use]
pub fn count_gc_case_insensitive(seq: &[u8]) -> u64 {
    let case = u8x16::splat(0x20);
    let g_lower = u8x16::splat(b'g');
    let c_lower = u8x16::splat(b'c');
    let mut count = 0u64;
    let chunks = seq.chunks_exact(16);
    let tail = chunks.remainder();
    for chunk in chunks {
        let v = u8x16::new(chunk.try_into().unwrap()) | case;
        // OR the two per-lane masks: no byte equals both `g` (0x67) and `c`
        // (0x63), so the masks are disjoint and a bitwise OR preserves the
        // total count without double-counting any lane.
        let gc_mask = v.simd_eq(g_lower).to_bitmask() | v.simd_eq(c_lower).to_bitmask();
        count += u64::from(gc_mask.count_ones());
    }
    for &b in tail {
        let folded = b | 0x20;
        if folded == b'g' || folded == b'c' {
            count += 1;
        }
    }
    count
}

/// Count `N` and `n` bases in `seq`.
///
/// Non-N bytes do not match; IUPAC ambiguity codes (`S`, `W`, `Y`, ...) are
/// treated as non-N, matching the behavior of the scalar
/// `match b { b'N' | b'n' => … }` form it replaces. Case-folding uses `| 0x20`
/// so both uppercase `N` and lowercase `n` are counted in a single compare.
#[must_use]
pub fn count_n_case_insensitive(seq: &[u8]) -> u64 {
    let case = u8x16::splat(0x20);
    let n_lower = u8x16::splat(b'n');
    let mut count = 0u64;
    let chunks = seq.chunks_exact(16);
    let tail = chunks.remainder();
    for chunk in chunks {
        let v = u8x16::new(chunk.try_into().unwrap()) | case;
        count += u64::from(v.simd_eq(n_lower).to_bitmask().count_ones());
    }
    for &b in tail {
        if b | 0x20 == b'n' {
            count += 1;
        }
    }
    count
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Scalar reference implementations ──────────────────────────────────────
    // Each SIMD kernel is tested by comparing its output against a trivial
    // byte-at-a-time reference on a variety of inputs.

    fn scalar_count_ge_q(qual: &[u8], threshold: u8) -> u64 {
        qual.iter().filter(|&&q| q >= threshold).count() as u64
    }

    fn scalar_count_lt_q(qual: &[u8], threshold: u8) -> u64 {
        qual.iter().filter(|&&q| q < threshold).count() as u64
    }

    fn scalar_count_gc(seq: &[u8]) -> u64 {
        seq.iter().filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c')).count() as u64
    }

    fn scalar_count_n(seq: &[u8]) -> u64 {
        seq.iter().filter(|&&b| matches!(b, b'N' | b'n')).count() as u64
    }

    // Tiny deterministic PRNG (xorshift64*) — avoids pulling in `rand` for tests.
    struct XorShift(u64);
    impl XorShift {
        fn new(seed: u64) -> Self {
            Self(seed.max(1))
        }
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            self.0 = x;
            x.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn fill_byte(&mut self, buf: &mut [u8], alphabet: &[u8]) {
            for b in buf {
                let r = self.next_u64();
                // Riker enforces a 64-bit or wider target in `lib.rs`, so the
                // cast is lossless. Allow clippy's pedantic truncation lint.
                #[allow(clippy::cast_possible_truncation)]
                let idx = (r as usize) % alphabet.len();
                *b = alphabet[idx];
            }
        }
    }

    // ── count_bases_ge_q / count_bases_lt_q ────────────────────────────────────

    #[test]
    fn count_bases_ge_q_empty() {
        assert_eq!(count_bases_ge_q(&[], 20), 0);
    }

    #[test]
    fn count_bases_ge_q_all_pass() {
        let qual = vec![40u8; 100];
        assert_eq!(count_bases_ge_q(&qual, 20), 100);
    }

    #[test]
    fn count_bases_ge_q_all_fail() {
        let qual = vec![10u8; 100];
        assert_eq!(count_bases_ge_q(&qual, 20), 0);
    }

    #[test]
    fn count_bases_ge_q_boundary_values() {
        // threshold is 20: 19 excluded, 20 included.
        let qual: Vec<u8> = (0..64).map(|i| if i % 2 == 0 { 19 } else { 20 }).collect();
        assert_eq!(count_bases_ge_q(&qual, 20), 32);
    }

    #[test]
    fn count_bases_ge_q_matches_scalar_on_chunk_boundaries() {
        let mut rng = XorShift::new(0xCAFE_F00D);
        // Every length across the SIMD chunk boundary, plus thresholds at 0
        // (degenerate "all pass") and 255 (degenerate "none pass") to pin
        // down the mask-extremes explicitly.
        for len in [0usize, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 150, 300] {
            let mut qual = vec![0u8; len];
            rng.fill_byte(&mut qual, &(0u8..=60).collect::<Vec<u8>>());
            for t in [0u8, 1, 10, 20, 30, 40, 60, 255] {
                assert_eq!(
                    count_bases_ge_q(&qual, t),
                    scalar_count_ge_q(&qual, t),
                    "len={len} t={t}",
                );
            }
        }
    }

    #[test]
    fn count_bases_lt_q_empty() {
        assert_eq!(count_bases_lt_q(&[], 20), 0);
    }

    #[test]
    fn count_bases_lt_q_matches_scalar_on_chunk_boundaries() {
        let mut rng = XorShift::new(0xDEAD_BEEF);
        for len in [0usize, 1, 15, 16, 17, 31, 32, 33, 150] {
            let mut qual = vec![0u8; len];
            rng.fill_byte(&mut qual, &(0u8..=60).collect::<Vec<u8>>());
            for t in [0u8, 1, 20, 30, 255] {
                assert_eq!(
                    count_bases_lt_q(&qual, t),
                    scalar_count_lt_q(&qual, t),
                    "len={len} t={t}",
                );
            }
        }
    }

    #[test]
    fn count_bases_ge_q_and_lt_q_partition_the_slice() {
        let mut rng = XorShift::new(0x1234_5678);
        let mut qual = vec![0u8; 200];
        rng.fill_byte(&mut qual, &(0u8..=60).collect::<Vec<u8>>());
        for t in [0u8, 1, 20, 30, 40] {
            assert_eq!(count_bases_ge_q(&qual, t) + count_bases_lt_q(&qual, t), qual.len() as u64,);
        }
    }

    // ── count_gc_case_insensitive ──────────────────────────────────────────────

    #[test]
    fn count_gc_empty() {
        assert_eq!(count_gc_case_insensitive(&[]), 0);
    }

    #[test]
    fn count_gc_case_insensitive_basic() {
        assert_eq!(count_gc_case_insensitive(b"ACGT"), 2);
        assert_eq!(count_gc_case_insensitive(b"acgt"), 2);
        assert_eq!(count_gc_case_insensitive(b"AaCcGgTt"), 4);
        assert_eq!(count_gc_case_insensitive(b"NNNNNN"), 0);
        assert_eq!(count_gc_case_insensitive(b""), 0);
    }

    #[test]
    fn count_gc_case_insensitive_rejects_ambiguity_codes() {
        // `S` encodes "G or C" in IUPAC but is not literally G/C; keep behavior
        // identical to the scalar `match b'G' | b'C' | b'g' | b'c'` it replaces.
        assert_eq!(count_gc_case_insensitive(b"SSSS"), 0);
        assert_eq!(count_gc_case_insensitive(b"KKKK"), 0);
    }

    #[test]
    fn count_gc_all_gc_boundary_lengths() {
        // Exercises the full-mask path (every lane of the u8x16 compare
        // set) at lengths that span the chunk boundary.
        for len in [15usize, 16, 17] {
            let seq: Vec<u8> = (0..len).map(|i| if i % 2 == 0 { b'G' } else { b'C' }).collect();
            assert_eq!(count_gc_case_insensitive(&seq), len as u64, "len={len}");
        }
        // Also verify the all-lowercase path, which exercises `| 0x20`
        // being a no-op.
        for len in [15usize, 16, 17] {
            let seq = vec![b'g'; len];
            assert_eq!(count_gc_case_insensitive(&seq), len as u64, "len={len} (all g)");
        }
    }

    #[test]
    fn count_gc_matches_scalar_on_chunk_boundaries() {
        let mut rng = XorShift::new(0xAA55_AA55);
        // Includes high-bit bytes (\x80, \xFF) to confirm the `| 0x20`
        // case-fold never produces a spurious match against `g`/`c`.
        let alphabet: &[u8] = b"ACGTNacgtnSKYMRW\x00\x80\xFF";
        for len in [0usize, 1, 15, 16, 17, 31, 32, 33, 64, 150, 300] {
            let mut seq = vec![0u8; len];
            rng.fill_byte(&mut seq, alphabet);
            assert_eq!(count_gc_case_insensitive(&seq), scalar_count_gc(&seq), "len={len}");
        }
    }

    // ── count_n_case_insensitive ───────────────────────────────────────────────

    #[test]
    fn count_n_empty() {
        assert_eq!(count_n_case_insensitive(&[]), 0);
    }

    #[test]
    fn count_n_case_insensitive_basic() {
        assert_eq!(count_n_case_insensitive(b"NNNN"), 4);
        assert_eq!(count_n_case_insensitive(b"nnnn"), 4);
        assert_eq!(count_n_case_insensitive(b"NnNn"), 4);
        assert_eq!(count_n_case_insensitive(b"ACGT"), 0);
    }

    #[test]
    fn count_n_matches_scalar_on_chunk_boundaries() {
        let mut rng = XorShift::new(0x3141_5926);
        // Same high-bit-byte guard as the GC test to confirm the `| 0x20`
        // fold never produces a spurious match against `n`.
        let alphabet: &[u8] = b"ACGTNacgtn\x00\x80\xFF";
        for len in [0usize, 1, 15, 16, 17, 31, 32, 33, 64, 150] {
            let mut seq = vec![0u8; len];
            rng.fill_byte(&mut seq, alphabet);
            assert_eq!(count_n_case_insensitive(&seq), scalar_count_n(&seq), "len={len}");
        }
    }
}
