//! Shared, test-only builders and assert helpers.
//!
//! This module is compiled only under `cfg(test)` (for the lib's own unit tests)
//! or the `test-support` feature (which `cargo ci-test` / `cargo ci-lint` turn on
//! for the separate integration-test crate). Nothing here ships in release builds.
//!
//! Its purpose is to give both test worlds *one* vocabulary for constructing test
//! inputs, replacing the record/FASTA/CIGAR builders that had drifted into a
//! dozen incompatible copies across `tests/` and inline `#[cfg(test)]` modules.
//! The design mirrors fgbio's `SamBuilder`/`ReferenceSetBuilder` ergonomics
//! (name only what differs) adapted to Rust via fluent setters over spec structs.

pub mod bed_builder;
pub mod cigar;
pub mod fasta_builder;
pub mod runner;
pub mod sam_builder;
pub mod vcf_builder;

pub use bed_builder::BedBuilder;
pub use cigar::ParsedCigar;
pub use fasta_builder::{DEFAULT_CONTIGS, FastaBuilder, default_contig_id};
pub use runner::drive_collector;
pub use sam_builder::{
    Adapters, IntoRecords, PairBuilder, ReadBuilder, SamBuilder, SortOrder, coord_builder, pair,
    read,
};
pub use vcf_builder::{InfoType, InfoValue, VariantEntry, VcfBuilder};

/// Assert two `f64` values are within `eps` (default `1e-6`) of each other.
///
/// The bound is inclusive — `|a - b| <= eps` — so a difference of exactly `eps` passes.
///
/// Usable from both worlds — unit tests via `crate::assert_close!`, integration
/// tests via `riker_lib::assert_close!` — replacing the hand-rolled
/// `(a - b).abs() < eps` (and the divergent `assert_float_eq!` / `assert_f64_eq!`)
/// scattered across the suite.
#[macro_export]
macro_rules! assert_close {
    ($a:expr, $b:expr $(,)?) => {
        $crate::assert_close!($a, $b, 1e-6_f64)
    };
    ($a:expr, $b:expr, $eps:expr $(,)?) => {{
        let a: f64 = $a;
        let b: f64 = $b;
        let eps: f64 = $eps;
        assert!(
            (a - b).abs() <= eps,
            "assert_close failed: |{a} - {b}| = {} > {eps}",
            (a - b).abs()
        );
    }};
}

#[cfg(test)]
mod tests {
    #[test]
    fn accepts_values_within_epsilon() {
        assert_close!(1.0, 1.0 + 1e-9);
        assert_close!(100.0, 100.4, 0.5);
    }

    #[test]
    #[should_panic(expected = "assert_close failed")]
    fn rejects_values_outside_epsilon() {
        assert_close!(1.0, 2.0);
    }
}
