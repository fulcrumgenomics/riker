#![deny(unsafe_code)]
#![allow(clippy::cast_precision_loss)]

use std::io::Write;
use std::time::Instant;

use anyhow::Result;
use chrono::Local;
use clap::Parser;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use env_logger::Env;
use riker_lib::commands::alignment::Alignment;
use riker_lib::commands::basic::Basic;
use riker_lib::commands::command::Command;
use riker_lib::commands::docs::Docs;
use riker_lib::commands::error::Error;
use riker_lib::commands::gcbias::GcBias;
use riker_lib::commands::hybcap::HybCap;
use riker_lib::commands::isize::InsertSize;
use riker_lib::commands::multi::Multi;
use riker_lib::commands::wgs::Wgs;

#[global_allocator]
static GLOBAL: rpmalloc::RpMalloc = rpmalloc::RpMalloc;

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default())
    .valid(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .invalid(AnsiColor::Yellow.on_default().effects(Effects::BOLD));

/// Fast sequencing QC metrics toolkit.
///
/// Riker collects alignment, insert-size, GC-bias, and whole-genome
/// coverage metrics from SAM, BAM, and CRAM files. Output is tab-separated with
/// lowercase snake_case headers and no metadata lines.
#[derive(Parser)]
#[command(
    name = "riker",
    version = env!("CARGO_PKG_VERSION"),
    long_about,
    styles = STYLES,
    after_long_help = "Run 'riker <command> --help' for detailed options on each command."
)]
struct Cli {
    /// Enable verbose logging.
    #[arg(long, global = true)]
    verbose: bool,

    #[command(subcommand)]
    command: Subcommand,
}

/// All riker subcommands.
#[derive(clap::Subcommand)]
enum Subcommand {
    Alignment(Alignment),
    Basic(Basic),
    Docs(Docs),
    Error(Error),
    Gcbias(GcBias),
    Hybcap(HybCap),
    Isize(InsertSize),
    Multi(Box<Multi>),
    Wgs(Wgs),
}

impl Command for Subcommand {
    fn execute(&self) -> Result<()> {
        match self {
            Subcommand::Alignment(c) => c.execute(),
            Subcommand::Basic(c) => c.execute(),
            Subcommand::Docs(c) => c.execute(),
            Subcommand::Error(c) => c.execute(),
            Subcommand::Gcbias(c) => c.execute(),
            Subcommand::Hybcap(c) => c.execute(),
            Subcommand::Isize(c) => c.execute(),
            Subcommand::Multi(c) => c.execute(),
            Subcommand::Wgs(c) => c.execute(),
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let log_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level))
        .format(|buf, record| {
            let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S");
            writeln!(buf, "[{timestamp} {} {}] {}", record.level(), record.target(), record.args())
        })
        .init();

    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    log::info!("Riker by Fulcrum Genomics - https://www.github.com/fulcrumgenomics/riker");
    log::info!("Executing: {cmdline}");

    let start = Instant::now();
    let result = cli.command.execute();
    let elapsed = start.elapsed();
    let minutes = elapsed.as_secs() / 60;
    let seconds = elapsed.as_secs() % 60;

    match &result {
        Ok(()) => log::info!("Successfully completed execution in {minutes}m:{seconds:02}s."),
        Err(_) => log::error!("Execution failed after {minutes}m:{seconds:02}s."),
    }

    result
}
