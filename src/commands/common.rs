use std::path::PathBuf;

pub const CORE_OPTIONS: &str = "Core Options";
pub const EXTENDED_OPTIONS: &str = "Extended Options";

/// Options for specifying the input alignment file.
#[derive(clap::Args, Debug, Clone)]
#[command()]
pub struct InputOptions {
    /// Input SAM, BAM, or CRAM file.
    #[arg(short = 'i', long, value_name = "INPUT")]
    pub input: PathBuf,
}

/// Options for specifying the output path prefix.
#[derive(clap::Args, Debug, Clone)]
#[command()]
pub struct OutputOptions {
    /// Output path prefix.
    #[arg(short = 'o', long, value_name = "PREFIX")]
    pub output: PathBuf,
}

/// Options for specifying a required reference FASTA.
#[derive(clap::Args, Debug, Clone)]
#[command()]
pub struct ReferenceOptions {
    /// Reference FASTA file (must be indexed with .fai). Required.
    #[arg(short = 'r', long, value_name = "FASTA")]
    pub reference: PathBuf,
}

/// Options for specifying an optional reference FASTA.
#[derive(clap::Args, Debug, Clone)]
#[command()]
pub struct OptionalReferenceOptions {
    /// Reference FASTA file (must be indexed with .fai). Optional unless CRAM input supplied.
    #[arg(short = 'r', long, value_name = "FASTA")]
    pub reference: Option<PathBuf>,
}
