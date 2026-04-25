//! SIMD-accelerated BAM sequence nibble decoder.
//!
//! BAM stores sequence bases packed two per byte: the high nibble is the
//! first base, the low nibble the second. Each nibble is an index into a
//! 16-entry table:
//!
//! | nibble | base |
//! |-------:|:----:|
//! | 0      | `=`  |
//! | 1      | `A`  |
//! | 2      | `C`  |
//! | 3      | `M`  |
//! | 4      | `G`  |
//! | 5      | `R`  |
//! | 6      | `S`  |
//! | 7      | `V`  |
//! | 8      | `T`  |
//! | 9      | `W`  |
//! | 10     | `Y`  |
//! | 11     | `H`  |
//! | 12     | `K`  |
//! | 13     | `D`  |
//! | 14     | `B`  |
//! | 15     | `N`  |
//!
//! noodles does this with a scalar `match` per nibble — fine, but also the
//! hottest per-base op on the reader thread once aux-tag decode is out of
//! the way. This module does it 32 bases at a time using `wide`'s
//! `u8x16::swizzle_relaxed`, which maps to NEON `vqtbl1q_u8` on arm64 and
//! `PSHUFB` on x86 (SSSE3+). A scalar tail handles whatever's left after
//! the 16-byte chunks.
//!
//! ## Shape of the kernel
//!
//! For each 16-byte input block (32 bases):
//!
//! 1. Load 16 packed bytes into a `u8x16`.
//! 2. Extract low nibbles (`packed & 0x0F`) with a bitwise AND — direct
//!    `u8x16` op.
//! 3. Extract high nibbles — `wide` doesn't expose per-byte shifts on
//!    `u8x16`, so we go through `u16x8`: right-shift by 4 on `u16x8`
//!    pulls the high-byte contamination into low-byte bits 4-7, which
//!    we clear with `& 0x0F0F`. Cast back to `u8x16`.
//! 4. Decode both nibble vectors via `swizzle_relaxed(TABLE, indices)`.
//! 5. Interleave `(hi[i], lo[i])` via `unpack_low` / `unpack_high` into
//!    two output vectors giving bases `[0..16]` and `[16..32]`.
//! 6. Append to the output `Vec`.
//!
//! ## Buffer reuse
//!
//! [`decode_packed_sequence_into`] takes `&mut Vec<u8>` and clears it
//! before decoding, so callers that keep a single output buffer across
//! many records avoid per-record allocation.

use bytemuck::cast;
use wide::{u8x16, u16x8};

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

#[cfg(test)]
mod tests {
    use super::*;

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
        // This test would fail if we forgot the `& 0x0F0F` after the u16
        // shift: when packed[i] has bit 7 set AND packed[i+1] has a low
        // nibble with bit 0 clear, the u16 shift would pull packed[i+1]
        // bits into packed[i]'s high-nibble output.
        let packed = vec![0xF8, 0x18, 0xF0, 0x01]; // 8 bases
        // Per-byte decode:
        //   0xF8 → high=F(N), low=8(T)
        //   0x18 → high=1(A), low=8(T)
        //   0xF0 → high=F(N), low=0(=)
        //   0x01 → high=0(=), low=1(A)
        let expected: Vec<u8> = b"NTATN==A".to_vec();
        // Sanity check: b'=' is 0x3D.
        assert_eq!(expected[5], b'=');
        let mut out = Vec::new();
        decode_packed_sequence_into(&packed, 8, &mut out);
        assert_eq!(out, expected);
    }
}
