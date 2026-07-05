//! Reusable SIMD kernels for riker's hot byte-level loops, built on the [`wide`] crate.
//!
//! Each kernel processes a fixed-width chunk of its input in SIMD and handles the
//! sub-chunk tail with a plain scalar loop. Most kernels are intentionally general
//! — "count bases ≥ Q over a slice", "GC count of a slice" — rather than tailored
//! to any one command, so a new call site can adopt them without growing the module.
//! A few are command-specific where the SIMD shape is inseparable from the caller's
//! per-position decision (e.g. [`pileup_depth_alone`], the `wgs` coverage kernel);
//! these live here to keep all `wide`-based vectorization in one module.
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

use bytemuck::cast;
use wide::{u8x16, u16x8};

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

/// Apply BQ-filtered, cap-saturated depth increments over a contiguous run of
/// pileup positions — the `AloneAction` fast path in `wgs::walk_depth`.
///
/// `qual` and `depth` are parallel slices of equal length. For each position `i`:
///
/// * `qual[i] < min_bq`         → a base-quality exclusion; `depth[i]` untouched.
/// * else `depth[i] >= cap`     → a capped exclusion; `depth[i]` untouched.
/// * else                       → `depth[i] += 1`.
///
/// Returns `(bases_excl_baseq, bases_excl_capped)`. The increment is applied only
/// where `depth[i] < cap`, so `depth[i] + 1 <= cap <= u16::MAX` and no overflow
/// is possible. Qualities are raw Phred values, not Phred+33.
///
/// This is the vectorized twin of the scalar per-position decision in
/// `wgs::WgsCollector::walk_depth`'s overlapping-mate branch; the two must stay
/// in lock-step (the byte-identical fixture run covers both).
#[must_use]
pub fn pileup_depth_alone(qual: &[u8], depth: &mut [u16], min_bq: u8, cap: u16) -> (u64, u64) {
    debug_assert_eq!(qual.len(), depth.len(), "qual and depth slices must be parallel");
    let n = qual.len().min(depth.len());

    let min_bq_v = u16x8::splat(u16::from(min_bq));
    let cap_v = u16x8::splat(cap);
    let one = u16x8::splat(1);
    let zero = u8x16::splat(0);

    let mut excl_baseq = 0u64;
    let mut excl_capped = 0u64;

    // 16 positions per iteration: quals load as u8x16 and widen into two u16x8
    // halves (interleave with a zero vector — little-endian puts the zero in the
    // high byte); depth is two native u16x8 vectors.
    let simd_len = n - (n % 16);
    let mut i = 0;
    while i < simd_len {
        let q16 = u8x16::new(qual[i..i + 16].try_into().unwrap());
        let q_lo: u16x8 = cast(u8x16::unpack_low(q16, zero));
        let q_hi: u16x8 = cast(u8x16::unpack_high(q16, zero));
        let d_lo = u16x8::new(depth[i..i + 8].try_into().unwrap());
        let d_hi = u16x8::new(depth[i + 8..i + 16].try_into().unwrap());

        let (d_lo2, bq_lo, cap_lo) = pileup_step(q_lo, d_lo, min_bq_v, cap_v, one);
        let (d_hi2, bq_hi, cap_hi) = pileup_step(q_hi, d_hi, min_bq_v, cap_v, one);

        let lo_arr: [u16; 8] = cast(d_lo2);
        let hi_arr: [u16; 8] = cast(d_hi2);
        depth[i..i + 8].copy_from_slice(&lo_arr);
        depth[i + 8..i + 16].copy_from_slice(&hi_arr);

        excl_baseq += bq_lo + bq_hi;
        excl_capped += cap_lo + cap_hi;
        i += 16;
    }

    // Scalar tail: the trailing 0..15 positions.
    for j in simd_len..n {
        let q = qual[j];
        if q < min_bq {
            excl_baseq += 1;
        } else if depth[j] < cap {
            depth[j] += 1;
        } else {
            excl_capped += 1;
        }
    }

    (excl_baseq, excl_capped)
}

/// One 8-lane pileup step: given widened quals `q` and current depths `d`,
/// return the updated depths plus the (BQ-fail, capped) lane counts.
#[inline]
fn pileup_step(q: u16x8, d: u16x8, min_bq_v: u16x8, cap_v: u16x8, one: u16x8) -> (u16x8, u64, u64) {
    // `wide`'s u16x8 exposes only simd_gt/eq, so express `x < y` as `y > x`.
    let bq_fail = min_bq_v.simd_gt(q); // min_bq > q ⟺ q < min_bq → BQ exclusion
    let below_cap = cap_v.simd_gt(d); // cap > depth ⟺ depth < cap → can take a read
    let pass = !bq_fail; // q >= min_bq
    let incr = pass & below_cap; // BQ-passing and below the cap → increment
    let capped = pass & !below_cap; // BQ-passing but already at the cap
    let d2 = d + (incr & one);
    let baseq = u64::from(bq_fail.to_bitmask().count_ones());
    let cap_excl = u64::from(capped.to_bitmask().count_ones());
    (d2, baseq, cap_excl)
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

// ─── BAM sequence nibble decoder ────────────────────────────────────────────────
//
// BAM stores sequence bases packed two per byte: the high nibble is the first base,
// the low nibble the second, each a 4-bit index into the 16-entry IUPAC table
// `NIBBLE_TO_BASE`. noodles decodes this with a scalar `match` per nibble; this does
// it 32 bases at a time via `wide`'s `u8x16::swizzle_relaxed` (NEON `vqtbl1q_u8` on
// arm64, `PSHUFB` on x86 SSSE3+), with a scalar tail for the leftover bytes. Per
// 16-byte block: mask the low nibbles directly; shift the high nibbles through
// `u16x8 >> 4` then `& 0x0F0F` (wide has no per-byte shift on `u8x16`); swizzle both
// through the table; interleave via unpack_low/high. `decode_packed_sequence_into`
// clears and reuses its `&mut Vec<u8>` so hot callers avoid per-record allocation.

/// BAM nibble → ASCII base lookup table. Indexing by `nibble & 0x0F`
/// yields the corresponding IUPAC code (case-insensitive).
pub(crate) const NIBBLE_TO_BASE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// Decode `base_count` bases from a `ceil(base_count / 2)`-byte packed
/// BAM sequence buffer into `dst`. The output is ASCII, one byte per base.
///
/// The destination `Vec` is cleared and re-populated. Existing capacity
/// is retained across calls.
///
/// # Panics
/// Panics in debug builds if `packed.len() < (base_count + 1) / 2`.
pub fn decode_packed_sequence_into(packed: &[u8], base_count: usize, dst: &mut Vec<u8>) {
    dst.clear();
    if base_count == 0 {
        return;
    }
    debug_assert!(
        packed.len() >= base_count.div_ceil(2),
        "packed buffer too short for {base_count} bases: have {}",
        packed.len()
    );

    dst.reserve(base_count);

    // SIMD path: 16 packed bytes → 32 bases per iteration.
    let table = u8x16::new(NIBBLE_TO_BASE);
    let mask_low_nibble = u8x16::splat(0x0F);
    // For the high-nibble step, we shift as u16x8 then mask per byte.
    let mask_both_low_nibbles_u16 = u16x8::splat(0x0F0F);

    let chunks = packed.chunks_exact(16);
    let tail = chunks.remainder();

    // Track how many bases we've written so we can stop precisely at
    // `base_count` (odd-length sequences don't use the final low nibble).
    let mut written: usize = 0;

    for chunk in chunks {
        // SAFETY of the unwrap: chunks_exact yields exactly 16 bytes.
        let packed_v = u8x16::new(chunk.try_into().unwrap());

        let low_nibbles = packed_v & mask_low_nibble;

        // High nibbles: shift via u16x8, mask with 0x0F0F, cast back.
        let as_u16: u16x8 = cast(packed_v);
        let shifted: u16x8 = as_u16 >> 4;
        let masked = shifted & mask_both_low_nibbles_u16;
        let high_nibbles: u8x16 = cast(masked);

        let decoded_hi = table.swizzle_relaxed(high_nibbles);
        let decoded_lo = table.swizzle_relaxed(low_nibbles);

        // Interleave so the output is [hi[0], lo[0], hi[1], lo[1], …].
        let first_16 = u8x16::unpack_low(decoded_hi, decoded_lo);
        let second_16 = u8x16::unpack_high(decoded_hi, decoded_lo);

        // How many of these 32 bases are still needed?
        let remaining = base_count - written;
        let first_16_arr: [u8; 16] = cast(first_16);
        let second_16_arr: [u8; 16] = cast(second_16);
        if remaining >= 32 {
            dst.extend_from_slice(&first_16_arr);
            dst.extend_from_slice(&second_16_arr);
            written += 32;
        } else {
            // Final chunk: take exactly `remaining` bases.
            let take_first = remaining.min(16);
            dst.extend_from_slice(&first_16_arr[..take_first]);
            let take_second = remaining - take_first;
            dst.extend_from_slice(&second_16_arr[..take_second]);
            // A fully-consumed stream with `base_count % 32 != 0` still has
            // a scalar tail ahead (from the un-chunked remainder), but it
            // won't contribute any additional bases — skip it.
            return;
        }
    }

    // Scalar tail: decode the 0..15 leftover packed bytes.
    scalar_tail(tail, base_count - written, dst);
}

/// Handle the leftover packed bytes after the SIMD loop. Small enough that
/// there's no benefit in SIMDing — and we also have to respect `base_count`
/// for odd-length sequences where the final low nibble is unused.
#[inline]
fn scalar_tail(packed_tail: &[u8], mut remaining: usize, dst: &mut Vec<u8>) {
    for &byte in packed_tail {
        if remaining == 0 {
            break;
        }
        let hi = NIBBLE_TO_BASE[(byte >> 4) as usize];
        dst.push(hi);
        remaining -= 1;
        if remaining == 0 {
            break;
        }
        let lo = NIBBLE_TO_BASE[(byte & 0x0F) as usize];
        dst.push(lo);
        remaining -= 1;
    }
}

/// Scalar reference implementation, kept for parity testing and as a
/// readable description of the intended semantics.
#[cfg(test)]
pub(crate) fn decode_packed_sequence_into_scalar(
    packed: &[u8],
    base_count: usize,
    dst: &mut Vec<u8>,
) {
    dst.clear();
    dst.reserve(base_count);
    scalar_tail(packed, base_count, dst);
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

    /// Scalar reference for `pileup_depth_alone`: returns the mutated depths and
    /// the (baseq-fail, capped) counts.
    fn scalar_pileup_alone(
        qual: &[u8],
        depth: &[u16],
        min_bq: u8,
        cap: u16,
    ) -> (Vec<u16>, u64, u64) {
        let mut d = depth.to_vec();
        let (mut bq, mut capped) = (0u64, 0u64);
        for (i, &q) in qual.iter().enumerate() {
            if q < min_bq {
                bq += 1;
            } else if d[i] < cap {
                d[i] += 1;
            } else {
                capped += 1;
            }
        }
        (d, bq, capped)
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

    // ── pileup_depth_alone ─────────────────────────────────────────────────────

    #[test]
    fn pileup_depth_alone_empty_is_noop() {
        let mut depth: [u16; 0] = [];
        assert_eq!(pileup_depth_alone(&[], &mut depth, 20, 250), (0, 0));
    }

    #[test]
    fn pileup_depth_alone_partitions_baseq_cap_and_increment() {
        // q<20 → baseq; q>=20 with depth<cap → +1; q>=20 with depth==cap → capped.
        let qual = [10u8, 30, 30, 40];
        let mut depth = [5u16, 7, 250, 0];
        let cap = 250;
        let (bq, capped) = pileup_depth_alone(&qual, &mut depth, 20, cap);
        assert_eq!(depth, [5, 8, 250, 1]); // pos0 baseq (untouched), pos2 at cap (untouched)
        assert_eq!(bq, 1);
        assert_eq!(capped, 1);
    }

    #[test]
    fn pileup_depth_alone_matches_scalar_across_boundaries() {
        let mut rng = XorShift::new(0xDA7A_5EED);
        // Lengths straddling the 16-wide SIMD block; depths span below/at/above
        // the cap so the cap and BQ mask lanes are all exercised.
        for len in [0usize, 1, 8, 15, 16, 17, 31, 32, 33, 100, 151, 300] {
            for _ in 0..16 {
                let mut qual = vec![0u8; len];
                rng.fill_byte(&mut qual, &(0u8..=60).collect::<Vec<u8>>());
                let depth: Vec<u16> =
                    (0..len).map(|_| u16::try_from(rng.next_u64() % 255).unwrap()).collect();
                for &min_bq in &[0u8, 20, 60] {
                    for &cap in &[1u16, 100, 250, 254] {
                        let (exp_d, exp_bq, exp_cap) =
                            scalar_pileup_alone(&qual, &depth, min_bq, cap);
                        let mut got_d = depth.clone();
                        let (got_bq, got_cap) = pileup_depth_alone(&qual, &mut got_d, min_bq, cap);
                        assert_eq!(got_d, exp_d, "depth len={len} min_bq={min_bq} cap={cap}");
                        assert_eq!(got_bq, exp_bq, "baseq len={len} min_bq={min_bq} cap={cap}");
                        assert_eq!(got_cap, exp_cap, "capped len={len} min_bq={min_bq} cap={cap}");
                    }
                }
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

    // ── BAM sequence nibble decoder ──

    /// Pack an ASCII base string into BAM 4-bit-packed bytes (inverse of
    /// the decoder). Used for test data construction.
    fn pack_bases(ascii: &[u8]) -> Vec<u8> {
        let mut packed = Vec::with_capacity(ascii.len().div_ceil(2));
        for chunk in ascii.chunks(2) {
            let hi = base_to_nibble(chunk[0]);
            let lo = if chunk.len() == 2 { base_to_nibble(chunk[1]) } else { 0 };
            packed.push((hi << 4) | lo);
        }
        packed
    }

    fn base_to_nibble(b: u8) -> u8 {
        u8::try_from(NIBBLE_TO_BASE.iter().position(|&x| x == b).unwrap())
            .expect("nibble index fits in u8 — table is 16 entries")
    }

    #[test]
    fn scalar_and_simd_agree_on_aligned_sequence() {
        let ascii: Vec<u8> = (0..256_usize).map(|i| NIBBLE_TO_BASE[i % 16]).collect();
        let packed = pack_bases(&ascii);

        let mut simd_out = Vec::new();
        decode_packed_sequence_into(&packed, ascii.len(), &mut simd_out);
        let mut scalar_out = Vec::new();
        decode_packed_sequence_into_scalar(&packed, ascii.len(), &mut scalar_out);

        assert_eq!(simd_out, ascii);
        assert_eq!(scalar_out, ascii);
    }

    #[test]
    fn scalar_and_simd_agree_on_odd_length() {
        // 33 bases → 17 packed bytes; the very last low-nibble is unused.
        let ascii: Vec<u8> = b"ACGTACGTNNNNACGTACGTNNNNACGTACGTA".to_vec();
        assert_eq!(ascii.len(), 33);
        let packed = pack_bases(&ascii);

        let mut simd_out = Vec::new();
        decode_packed_sequence_into(&packed, ascii.len(), &mut simd_out);
        assert_eq!(simd_out, ascii);
    }

    #[test]
    fn short_sequences_under_simd_width() {
        // < 16 packed bytes means the SIMD loop runs zero times; the
        // scalar tail has to produce the whole sequence.
        for len in [1usize, 2, 15, 16, 17, 31, 32, 33] {
            let ascii: Vec<u8> = (0..len).map(|i| NIBBLE_TO_BASE[i % 16]).collect();
            let packed = pack_bases(&ascii);
            let mut out = Vec::new();
            decode_packed_sequence_into(&packed, len, &mut out);
            assert_eq!(out, ascii, "mismatch at len={len}");
        }
    }

    #[test]
    fn empty_sequence_produces_empty_output() {
        let mut out = vec![1, 2, 3];
        decode_packed_sequence_into(&[], 0, &mut out);
        assert!(out.is_empty(), "empty input should clear the destination");
    }

    #[test]
    fn existing_destination_capacity_is_reused() {
        let ascii: Vec<u8> = (0..150).map(|i| NIBBLE_TO_BASE[i % 16]).collect();
        let packed = pack_bases(&ascii);
        let mut out: Vec<u8> = Vec::with_capacity(200);
        let cap_before = out.capacity();
        decode_packed_sequence_into(&packed, ascii.len(), &mut out);
        let cap_after = out.capacity();
        assert_eq!(out, ascii);
        // Reuse: we pre-sized to >= needed, so no reallocation.
        assert_eq!(cap_after, cap_before);
    }

    #[test]
    fn realistic_150bp_read() {
        // Canonical Illumina-ish length; checks the common hot-path case.
        let ascii: Vec<u8> = (0..150).map(|i| b"ACGT"[i % 4]).collect();
        let packed = pack_bases(&ascii);
        let mut out = Vec::new();
        decode_packed_sequence_into(&packed, ascii.len(), &mut out);
        assert_eq!(out, ascii);
    }

    /// Quick head-to-head timing of the SIMD decoder vs the scalar
    /// reference on a realistic Illumina-sized read.
    ///
    /// Not a real benchmark — just enough instrumentation to show the
    /// relative speedup during development. Run with:
    ///
    /// ```text
    /// cargo test --release -- --ignored --nocapture bench_simd_decode
    /// ```
    #[test]
    #[ignore = "perf instrumentation; run with --release --ignored"]
    fn bench_simd_decode() {
        use std::time::Instant;

        // ~Illumina read length: 150bp. 75 packed bytes per record.
        let ascii: Vec<u8> = (0..150).map(|i| b"ACGT"[i % 4]).collect();
        let packed = pack_bases(&ascii);
        let iters: u32 = 5_000_000;
        let mut out = Vec::with_capacity(ascii.len());

        // SIMD
        let start = Instant::now();
        for _ in 0..iters {
            decode_packed_sequence_into(&packed, ascii.len(), &mut out);
            std::hint::black_box(&out);
        }
        let simd = start.elapsed();

        // Scalar
        let start = Instant::now();
        for _ in 0..iters {
            decode_packed_sequence_into_scalar(&packed, ascii.len(), &mut out);
            std::hint::black_box(&out);
        }
        let scalar = start.elapsed();

        let simd_ns_per = simd.as_nanos() as f64 / f64::from(iters);
        let scalar_ns_per = scalar.as_nanos() as f64 / f64::from(iters);
        eprintln!(
            "150bp decode × {iters}:  simd {simd:?} ({simd_ns_per:.1} ns/rec)  scalar {scalar:?} ({scalar_ns_per:.1} ns/rec)  speedup {:.2}x",
            scalar_ns_per / simd_ns_per,
        );
    }

    #[test]
    fn high_nibble_mask_prevents_cross_byte_leak() {
        // This test would fail if we forgot the `& 0x0F0F` after the u16 shift: when
        // packed[i] has bit 7 set AND packed[i+1] has a low nibble with bit 0 clear, the
        // u16 shift would pull packed[i+1]'s bits into packed[i]'s high-nibble output.
        //
        // Pad to a full 16-byte chunk so the leak-prone bytes go through the SIMD path
        // (u16 shift + mask), not the scalar tail.
        let mut packed = vec![0xF8, 0x18, 0xF0, 0x01];
        packed.resize(16, 0x00); // padding bytes → nibble 0 → '=' bases
        // First four bytes decode as:
        //   0xF8 → high=F(N), low=8(T)
        //   0x18 → high=1(A), low=8(T)
        //   0xF0 → high=F(N), low=0(=)
        //   0x01 → high=0(=), low=1(A)
        let mut expected = b"NTATN==A".to_vec();
        expected.resize(32, b'='); // 12 padding bytes → 24 '=' bases
        // Sanity check: b'=' is 0x3D.
        assert_eq!(expected[5], b'=');
        let mut out = Vec::new();
        decode_packed_sequence_into(&packed, 32, &mut out);
        assert_eq!(out, expected);
    }
}
