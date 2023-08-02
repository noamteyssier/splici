use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Command,
}

#[derive(Subcommand)]
pub enum Command {
    /// Select Introns and Exons from a GTF and parse sequences
    /// from a fasta file
    Introns {
        /// Input GTF file to convert
        #[clap(short, long)]
        input: String,

        /// Input fasta file to query intervals from
        #[clap(short, long)]
        fasta: String,

        /// Output file to write to (default=stdout)
        #[clap(short, long)]
        output: Option<String>,

        /// Read length to extend introns by
        #[clap(short, long, default_value = "150")]
        extension: Option<usize>,
    },
}
