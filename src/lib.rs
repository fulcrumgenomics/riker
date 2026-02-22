#![deny(unsafe_code)]
#![allow(clippy::cast_precision_loss)]

#[cfg(any(target_pointer_width = "16", target_pointer_width = "32"))]
compile_error!("riker requires a 64-bit or wider platform");

pub mod collector;
pub mod commands;
pub mod counter;
pub mod fasta;
pub mod intervals;
pub mod math;
pub mod metrics;
pub mod overlapper;
pub mod plotting;
pub mod progress;
pub mod sam;
pub mod sequence_dict;
pub mod vcf;
pub mod version;
