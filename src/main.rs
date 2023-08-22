mod cli;
mod functions;
mod io;
mod types;
mod utils;
use std::time::SystemTime;

use anyhow::Result;
use bedrs::{types::StrandedGenomicIntervalSet, StrandedGenomicInterval};
use clap::Parser;
use cli::{Cli, Command};
use hashbrown::HashMap;

pub type Giv = StrandedGenomicInterval<usize>;
pub type GivSet = StrandedGenomicIntervalSet<usize>;
pub type NameMap = HashMap<Vec<u8>, usize>;
pub type IdMap = HashMap<usize, Vec<u8>>;

fn setup_logger() -> Result<()> {
    fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "[{} {}] {}",
                humantime::format_rfc3339_seconds(SystemTime::now()),
                record.level(),
                message
            ))
        })
        .level(log::LevelFilter::Debug)
        .chain(std::io::stderr())
        .apply()?;
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    setup_logger()?;
    match cli.command {
        Command::Introns {
            gtf,
            fasta,
            output,
            t2g,
            extension,
            num_threads,
            compression_level,
        } => {
            functions::run_introns(
                &gtf,
                &fasta,
                t2g,
                output,
                extension,
                num_threads,
                compression_level,
            )?;
        }
    }
    Ok(())
}
