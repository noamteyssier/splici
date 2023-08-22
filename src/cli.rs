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
        gtf: String,

        /// Input fasta file to query intervals from
        #[clap(short, long)]
        fasta: String,

        /// Output file to write to (default=stdout)
        #[clap(short, long)]
        output: Option<String>,

        /// Output Transcript to Gene mapping file
        #[clap(short, long, default_value = "t2g.tsv")]
        t2g: String,

        /// Read length to extend introns by
        #[clap(short, long, default_value = "150")]
        extension: Option<usize>,

        #[clap(short = 'j', long, value_parser)]
        /// Number of threads to use in gzip compression
        num_threads: Option<usize>,

        #[clap(short = 'Z', long, value_parser)]
        /// gzip compression level
        compression_level: Option<usize>,
    },
}
