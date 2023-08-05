use anyhow::Result;
use flate2::read::MultiGzDecoder;
use gzp::{
    deflate::Gzip,
    par::compress::{ParCompress, ParCompressBuilder},
    Compression,
};
use std::{
    fs::File,
    io::{stdout, BufRead, BufReader, BufWriter, Write},
};

/// Matches the output to a writer stream
pub fn match_output_stream(
    output: Option<String>,
    num_threads: Option<usize>,
    compression_level: Option<usize>,
) -> Result<Box<dyn Write>> {
    match output {
        Some(path) => {
            let file = File::create(&path)?;
            if path.ends_with(".gz") {
                let writer: ParCompress<Gzip> = ParCompressBuilder::new()
                    .num_threads(num_threads.unwrap_or(1))?
                    .compression_level(if let Some(level) = compression_level {
                        Compression::new(level as u32)
                    } else {
                        Compression::default()
                    })
                    .from_writer(file);
                let buffer = BufWriter::new(writer);
                Ok(Box::new(buffer))
            } else {
                let buffer = BufWriter::new(file);
                Ok(Box::new(buffer))
            }
        }
        None => {
            let buffer = BufWriter::new(stdout());
            Ok(Box::new(buffer))
        }
    }
}

pub fn match_input_gtf(input: &str) -> Result<Box<dyn BufRead>> {
    let handle = File::open(input)?;
    if input.ends_with(".gz") {
        let decoder = MultiGzDecoder::new(handle);
        let buffer = BufReader::new(decoder);
        Ok(Box::new(buffer))
    } else {
        let buffer = BufReader::new(handle);
        Ok(Box::new(buffer))
    }
}
