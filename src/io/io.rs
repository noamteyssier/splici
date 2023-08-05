use anyhow::Result;
use gzp::{
    deflate::Gzip,
    par::compress::{ParCompress, ParCompressBuilder},
    Compression,
};
use std::{
    fs::File,
    io::{stdout, Write},
};

/// Matches the output to a writer stream
pub fn match_output_stream(
    output: Option<String>,
    num_threads: Option<usize>,
    compression_level: Option<usize>,
) -> Result<Box<dyn Write>> {
    match output {
        Some(path) => {
            if path.ends_with(".gz") {
                let file = File::create(path)?;
                let writer: ParCompress<Gzip> = ParCompressBuilder::new()
                    .num_threads(num_threads.unwrap_or(1))?
                    .compression_level(if let Some(level) = compression_level {
                        Compression::new(level as u32)
                    } else {
                        Compression::default()
                    })
                    .from_writer(file);
                Ok(Box::new(writer))
            } else {
                Ok(Box::new(File::create(path)?))
            }
        }
        None => Ok(Box::new(stdout())),
    }
}
