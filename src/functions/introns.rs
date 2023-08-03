use crate::{
    types::ExonRecord,
    utils::{build_interval_set, flip_map, get_gene, get_introns, interval_to_region, build_exon_set, merge_interval_set, parse_exons},
};
use anyhow::Result;
use bedrs::{Container, GenomicInterval, GenomicIntervalSet};
use gtftools::GtfReader;
use hashbrown::HashMap;
use noodles::fasta::{fai, IndexedReader, io::BufReadSeek};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    str::from_utf8,
};

struct Splici {
    /// Stores a mapping of genome_id to genome_name
    genome_map: HashMap<usize, Vec<u8>>,

    /// Stores a mapping of genome_name to genome_id
    genome_names: HashMap<Vec<u8>, usize>,

    /// Stores a mapping of gene_id to gene_name
    gene_map: HashMap<usize, Vec<u8>>,

    /// Stores a mapping of gene_name to gene_id
    gene_names: HashMap<Vec<u8>, usize>,

    /// Stores a mapping of transcript_id to transcript_name
    transcript_map: HashMap<usize, Vec<u8>>,

    /// Stores a mapping of transcript_name to transcript_id
    transcript_names: HashMap<Vec<u8>, usize>,

    /// Stores a mapping of transcript_id to exon records
    transcript_records: HashMap<usize, Vec<ExonRecord>>,

    /// Stores a mapping of gene_id to a vector of intron records
    gene_introns: HashMap<usize, Vec<GenomicInterval<usize>>>,

    /// Stores a mapping of gene_id to merged intron sets
    merged_introns: HashMap<usize, GenomicIntervalSet<usize>>,
}
impl Splici {
    pub fn new() -> Self {
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
        }
    }

    pub fn parse_exons<R>(&mut self, reader: &mut GtfReader<R>) -> Result<()>
    where
        R: BufRead,
    {
        self.transcript_records = parse_exons(
            reader,
            &mut self.genome_names,
            &mut self.gene_names,
            &mut self.transcript_names,
        )?;
        self.flip_maps();
        Ok(())
    }

    fn flip_maps(&mut self) {
        self.genome_map = flip_map(&self.genome_names);
        self.gene_map = flip_map(&self.gene_names);
        self.transcript_map = flip_map(&self.transcript_names);
    }

    pub fn parse_introns(&mut self) {
        self.gene_introns = self
            .transcript_records
            .values()
            .map(|exon_vec| (get_gene(exon_vec), exon_vec))
            .map(|(gene_id, exon_vec)| (gene_id, build_exon_set(exon_vec)))
            .map(|(gene_id, exon_set)| (gene_id, get_introns(exon_set)))
            .fold(HashMap::new(), |mut gene_introns, (gene_id, intron_set)| {
                intron_set.iter().for_each(|intron| {
                    gene_introns
                        .entry(gene_id)
                        .or_insert_with(Vec::new)
                        .push(*intron);
                });
                gene_introns
            });
    }

    pub fn merge_introns(&mut self) {
        self.merged_introns = self.gene_introns
            .iter()
            .map(|(gene_id, intron_vec)| {
                (gene_id, build_interval_set(&intron_vec))
            })
            .map(|(gene_id, intron_set)| {
                (*gene_id, merge_interval_set(intron_set))
            })
            .collect();
    }

    fn write_introns<R: BufReadSeek>(&self, fasta: &mut IndexedReader<R>) {
        self.merged_introns
            .keys()
            .for_each(|gene_id| self.write_introns_for(*gene_id, fasta))
    }

    fn write_exons<R: BufReadSeek>(&self, fasta: &mut IndexedReader<R>) {
        self.transcript_records
            .keys()
            .for_each(|tx| self.write_exons_for(*tx, fasta))
    }

    fn write_introns_for<R: BufReadSeek>(&self, gene_id: usize, fasta: &mut IndexedReader<R>) {
        let intron_set = self.merged_introns.get(&gene_id).unwrap();
        let gene_name = self.gene_map.get(&gene_id).map(|x| from_utf8(x).unwrap()).unwrap();
        intron_set
            .records()
            .iter()
            .enumerate()
            .for_each(|(idx, x)| {
                let region = interval_to_region(x, &self.genome_map).unwrap();
                let query = fasta.query(&region).unwrap();
                let seq = from_utf8(query.sequence().as_ref()).unwrap();
                print!(">{gene_name}-I.{idx}\n{seq}\n");
            })
    }

    fn write_exons_for<R: BufReadSeek>(&self, tx: usize, fasta: &mut IndexedReader<R>) {
        let transcript_name = self.get_transcript_name(tx).expect("Could not get transcript name");
        let exon_set = self.get_exon_set(tx);
        let transcript_seq = self.build_transcript_sequence(exon_set, fasta);
        println!(">{}", transcript_name);
        println!("{}", from_utf8(&transcript_seq).unwrap());
    }

    fn get_transcript_name(&self, tx: usize) -> Option<&str> {
        self.transcript_map
            .get(&tx)
            .map(|x| from_utf8(x).unwrap())
    }

    fn get_exon_set(&self, tx: usize) -> GenomicIntervalSet<usize> {
        let mut exon_set: GenomicIntervalSet<usize> = self.transcript_records
            .get(&tx)
            .unwrap()
            .iter()
            .map(|x| x.into())
            .collect();
        exon_set.sort();
        exon_set
    }

    fn build_transcript_sequence<R: BufReadSeek>(&self, exon_set: GenomicIntervalSet<usize>, fasta: &mut IndexedReader<R>) -> Vec<u8> {
        let mut transcript_seq: Vec<u8> = Vec::new();
        exon_set.records().iter().for_each(|exon| {
            let exon_seq = self.get_exon_sequence(exon, fasta);
            transcript_seq.extend(&exon_seq);
        });
        transcript_seq
    }

    fn get_exon_sequence<R: BufReadSeek>(&self, exon: &GenomicInterval<usize>, fasta: &mut IndexedReader<R>) -> Vec<u8>  {
        let region = interval_to_region(exon, &self.genome_map).expect("Could not get exon sequence");
        let query = fasta.query(&region).expect("Could not query fasta");
        query.sequence().as_ref().to_vec()
    }

}

pub fn run_introns(
    gtf_path: &str,
    fasta_path: &str,
    output_path: Option<String>,
    extension: Option<usize>,
) -> Result<()> {
    let index_file_path = format!("{fasta_path}.fai");
    let index = fai::read(index_file_path)?;
    let reader = File::open(fasta_path).map(BufReader::new)?;
    let mut fasta = IndexedReader::new(reader, index);

    let handle = File::open(gtf_path).map(BufReader::new)?;
    let mut reader = GtfReader::from_bufread(handle);

    let mut splici = Splici::new();
    splici.parse_exons(&mut reader)?;
    splici.parse_introns();
    splici.merge_introns();
    splici.write_exons(&mut fasta);
    splici.write_introns(&mut fasta);

    Ok(())
}
