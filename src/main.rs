mod cli;
mod functions;
mod types;
mod utils;
use anyhow::Result;
use clap::Parser;
use cli::{Cli, Command};

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Introns {
            input,
            fasta,
            output,
            extension,
        } => {
            functions::run(&input, &fasta, output, extension)?;
        }
    }
    Ok(())
}
