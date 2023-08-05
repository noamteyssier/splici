use crate::{
    io::{match_input_gtf, match_output_stream},
    types::ExonRecord,
    utils::{
        build_exon_set, build_interval_set, flip_map, get_gene, get_introns, interval_to_region,
        merge_interval_set, parse_exons, reverse_complement,
    },
    Giv, GivSet, IdMap, NameMap,
};
use anyhow::Result;
use bedrs::{Container, Strand};
use gtftools::GtfReader;
use hashbrown::HashMap;
use log::info;
use noodles::fasta::{fai, io::BufReadSeek, IndexedReader};
use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
    str::from_utf8,
};

struct Splici {
    /// Stores a mapping of genome_id to genome_name
    genome_map: IdMap,

    /// Stores a mapping of gene_id to gene_name
    gene_map: IdMap,

    /// Stores a mapping of transcript_id to transcript_name
    transcript_map: IdMap,

    /// Stores a mapping of genome_name to genome_id
    genome_names: NameMap,

    /// Stores a mapping of gene_name to gene_id
    gene_names: NameMap,

    /// Stores a mapping of transcript_name to transcript_id
    transcript_names: NameMap,

    /// Stores a mapping of transcript_id to exon records
    transcript_records: HashMap<usize, Vec<ExonRecord>>,

    /// Stores a mapping of gene_id to a vector of intron records
    gene_introns: HashMap<usize, Vec<Giv>>,

    /// Stores a mapping of gene_id to merged intron sets
    merged_introns: HashMap<usize, GivSet>,

    /// The extension length for intronic regions
    extension: Option<usize>,
}
impl Splici {
    pub fn new(extension: Option<usize>) -> Self {
        Self {
            genome_map: HashMap::new(),
            gene_map: HashMap::new(),
            transcript_map: HashMap::new(),
            genome_names: HashMap::new(),
            gene_names: HashMap::new(),
            transcript_names: HashMap::new(),
            transcript_records: HashMap::new(),
            gene_introns: HashMap::new(),
            merged_introns: HashMap::new(),
            extension,
        }
    }

    /// Parses the exons from a GTF file into a `HashMap` of `transcript_id` to a Vec of `ExonRecords`
    pub fn parse_exons<R>(&mut self, reader: &mut GtfReader<R>) -> Result<()>
    where
        R: BufRead,
    {
        info!("Parsing exons from GTF");
        self.transcript_records = parse_exons(
            reader,
            &mut self.genome_names,
            &mut self.gene_names,
            &mut self.transcript_names,
        )?;

        if self.transcript_records.is_empty() {
            anyhow::bail!("No exons found in GTF (is it a properly formatted GTF?)");
        }

        self.flip_maps();
        let num_exons = self
            .transcript_records
            .values()
            .map(|exon_vec| exon_vec.len())
            .sum::<usize>();

        info!("Parsed {} genes", self.gene_names.len());
        info!("Parsed {} transcripts", self.transcript_names.len());
        info!("Parsed {} exons", num_exons);
        Ok(())
    }

    /// Flips the genome, gene, and transcript maps
    /// This is necessary because we need O(1) lookup in the opposite direction
    /// i.e. `genome_id` -> `genome_name` instead of `genome_name` -> `genome_id`
    fn flip_maps(&mut self) {
        self.genome_map = flip_map(&self.genome_names);
        self.gene_map = flip_map(&self.gene_names);
        self.transcript_map = flip_map(&self.transcript_names);
    }

    /// Parses the introns set for each transcript
    pub fn parse_introns(&mut self) -> Result<()> {
        let mut intron_count = 0;
        info!("Calculating intronic regions");
        for exon_vec in self.transcript_records.values() {
            let gene_id = get_gene(exon_vec);
            let exon_set = build_exon_set(exon_vec);
            let intron_set = get_introns(exon_set, self.extension)?;
            for intron in intron_set {
                self.gene_introns
                    .entry(gene_id)
                    .or_insert_with(Vec::new)
                    .push(intron);
                intron_count += 1;
            }
        }
        info!("Identified {} introns", intron_count);
        Ok(())
    }

    /// Merges the overlapping intronic regions for each gene
    pub fn merge_introns(&mut self) {
        info!("Merging overlapping intronic regions");
        let mut merged_count = 0;
        self.merged_introns = self
            .gene_introns
            .iter()
            .map(|(gene_id, intron_vec)| (gene_id, build_interval_set(intron_vec)))
            .map(|(gene_id, intron_set)| (*gene_id, merge_interval_set(intron_set)))
            .map(|(gene_id, intron_set)| {
                merged_count += intron_set.len();
                (gene_id, intron_set)
            })
            .collect();
        info!("Identified {} introns", merged_count);
    }

    /// Writes all intronic region sequences to stdout
    fn write_introns<R, W>(&self, fasta: &mut IndexedReader<R>, writer: &mut W) -> Result<()>
    where
        R: BufReadSeek,
        W: Write,
    {
        info!("Writing intronic regions");
        for gene_id in self.merged_introns.keys() {
            self.write_introns_for(*gene_id, fasta, writer)?;
        }
        info!("Finished writing intronic regions");
        Ok(())
    }

    /// Writes all concatenated exon transcripts to stdout
    fn write_exons<R, W>(&self, fasta: &mut IndexedReader<R>, writer: &mut W) -> Result<()>
    where
        R: BufReadSeek,
        W: Write,
    {
        info!("Writing exon transcripts");
        for transcript_id in self.transcript_records.keys() {
            self.write_exons_for(*transcript_id, fasta, writer)?;
        }
        info!("Finished writing exon transcripts");
        Ok(())
    }

    /// Writes the specific intronic region sequences for a given gene to stdout
    fn write_introns_for<R, W>(
        &self,
        gene_id: usize,
        fasta: &mut IndexedReader<R>,
        writer: &mut W,
    ) -> Result<()>
    where
        R: BufReadSeek,
        W: Write,
    {
        let intron_set = self.merged_introns.get(&gene_id).unwrap();
        let gene_name = from_utf8(self.gene_map.get(&gene_id).unwrap())?;
        let intron_iter = intron_set.records().iter().enumerate();
        for (idx, x) in intron_iter {
            let region = interval_to_region(x, &self.genome_map).unwrap();
            let query = fasta.query(&region).unwrap();
            let utf_seq = query.sequence().as_ref();
            let seq = match x.strand() {
                Strand::Reverse => reverse_complement(utf_seq),
                _ => utf_seq.to_vec(),
            };
            let seq = from_utf8(&seq).unwrap();
            write!(writer, ">{gene_name}-I.{idx}\n{seq}\n")?;
        }
        Ok(())
    }

    /// Writes the specific exon transcripts for a given transcript to stdout
    fn write_exons_for<R, W>(
        &self,
        tx: usize,
        fasta: &mut IndexedReader<R>,
        writer: &mut W,
    ) -> Result<()>
    where
        R: BufReadSeek,
        W: Write,
    {
        let transcript_name = self
            .get_transcript_name(tx)
            .expect("Could not get transcript name");
        let exon_set = self.get_exon_set(tx);
        let strand = exon_set.records()[0].strand();
        let transcript_seq = self.build_transcript_sequence(exon_set, fasta);
        let transcript_seq = match strand {
            Strand::Reverse => reverse_complement(&transcript_seq),
            _ => transcript_seq,
        };
        let transcript_str = from_utf8(&transcript_seq)?;
        write!(writer, ">{transcript_name}\n{transcript_str}\n")?;
        Ok(())
    }

    /// Convenience function to get the transcript name from the transcript id
    fn get_transcript_name(&self, tx: usize) -> Option<&str> {
        self.transcript_map.get(&tx).map(|x| from_utf8(x).unwrap())
    }

    /// Convenience function to get the exon set from the transcript id
    fn get_exon_set(&self, tx: usize) -> GivSet {
        let mut exon_set: GivSet = self
            .transcript_records
            .get(&tx)
            .unwrap()
            .iter()
            .map(std::convert::Into::into)
            .collect();
        exon_set.sort();
        exon_set
    }

    /// Builds the transcript sequence from the exon set
    fn build_transcript_sequence<R: BufReadSeek>(
        &self,
        exon_set: GivSet,
        fasta: &mut IndexedReader<R>,
    ) -> Vec<u8> {
        let mut transcript_seq: Vec<u8> = Vec::new();
        exon_set.records().iter().for_each(|exon| {
            let region = interval_to_region(exon, &self.genome_map).unwrap();
            let query = fasta.query(&region).unwrap();
            let seq = query.sequence().as_ref();
            transcript_seq.extend(seq);
        });
        transcript_seq
    }
}

pub fn run_introns(
    gtf_path: &str,
    fasta_path: &str,
    output_path: Option<String>,
    extension: Option<usize>,
    num_threads: Option<usize>,
    compression_level: Option<usize>,
) -> Result<()> {
    let index_file_path = format!("{fasta_path}.fai");
    let index = fai::read(index_file_path)?;
    let reader = File::open(fasta_path).map(BufReader::new)?;
    let mut fasta = IndexedReader::new(reader, index);

    // let handle = File::open(gtf_path).map(BufReader::new)?;
    let gtf_handle = match_input_gtf(gtf_path)?;
    let mut reader = GtfReader::from_bufread(gtf_handle);
    let mut output = match_output_stream(output_path, num_threads, compression_level)?;

    let mut splici = Splici::new(extension);
    splici.parse_exons(&mut reader)?;
    splici.parse_introns()?;
    splici.merge_introns();
    splici.write_exons(&mut fasta, &mut output)?;
    splici.write_introns(&mut fasta, &mut output)?;

    Ok(())
}
