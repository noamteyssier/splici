mod cli;
mod functions;
mod types;
mod utils;
use anyhow::Result;
use bedrs::{GenomicInterval, GenomicIntervalSet};
use clap::Parser;
use cli::{Cli, Command};
use hashbrown::HashMap;

pub type Giv = GenomicInterval<usize>;
pub type GivSet = GenomicIntervalSet<usize>;
pub type NameMap = HashMap<Vec<u8>, usize>;
pub type IdMap = HashMap<usize, Vec<u8>>;

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Introns {
            input,
            fasta,
            output,
            extension,
        } => {
            functions::run_introns(&input, &fasta, output, extension)?;
        }
    }
    Ok(())
}
