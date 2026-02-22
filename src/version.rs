#![allow(clippy::doc_markdown)] // Generated file contains OPT_LEVEL without backticks

use std::sync::LazyLock;

include!(concat!(env!("OUT_DIR"), "/built.rs"));

/// Version string from Cargo package metadata.
pub static VERSION: LazyLock<String> = LazyLock::new(|| PKG_VERSION.to_string());
