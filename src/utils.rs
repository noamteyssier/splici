use anyhow::{bail, Result};
use bedrs::{Container, Coordinates, GenomicInterval, GenomicIntervalSet, Internal, Merge};
use gtftools::GtfReader;
use hashbrown::HashMap;
use noodles::core::{Position, Region};
use std::{hash::Hash, io::BufRead, str::from_utf8};

use crate::types::ExonRecord;

/// Flips the K, V pairs in a HashMap
/// Used to reverse the name -> idx mappings to idx -> name
pub fn flip_map<K, V>(map: &HashMap<K, V>) -> HashMap<V, K>
where
    K: Eq + Hash + Clone,
    V: Eq + Hash + Clone,
{
    let mut flipped = HashMap::new();
    for (key, value) in map.iter() {
        flipped.insert(value.clone(), key.clone());
    }
    flipped
}

/// Converts a `GenomicInterval` to a `Region`
pub fn interval_to_region(
    giv: &GenomicInterval<usize>,
    genome_name_map: &HashMap<usize, Vec<u8>>,
) -> Result<Region> {
    let name = if let Some(name) = genome_name_map.get(&giv.chr()) {
        from_utf8(name)?.to_string()
    } else {
        bail!("Could not find genome name for id {}", giv.chr())
    };
    let start = Position::try_from(giv.start())?;
    let end = Position::try_from(giv.end())?;
    let region = Region::new(name, start..=end);
    Ok(region)
}

pub fn get_gene(exon_set: &[ExonRecord]) -> usize {
    exon_set[0].gene()
}

pub fn build_exon_set(exon_set: &[ExonRecord]) -> GenomicIntervalSet<usize> {
    let mut exon_set = exon_set
        .iter()
        .map(|exon| exon.into())
        .collect::<GenomicIntervalSet<usize>>();
    exon_set.sort();
    exon_set
}

pub fn build_interval_set(vec: &[GenomicInterval<usize>]) -> GenomicIntervalSet<usize> {
    let mut set = GenomicIntervalSet::from_iter(vec.iter().cloned());
    set.sort();
    set
}

pub fn get_introns(giv_set: GenomicIntervalSet<usize>) -> Vec<GenomicInterval<usize>> {
    giv_set
        .internal()
        .expect("could not parse introns")
        .collect()
}

pub fn merge_interval_set(giv_set: GenomicIntervalSet<usize>) -> GenomicIntervalSet<usize> {
    giv_set
        .merge()
        .expect("Could not merge interval set")
        .intervals()
        .into_iter()
        .cloned()
        .collect()
}

/// Parse the exons from the GTF file and return a HashMap of transcript id to
/// a Vec of `ExonRecords`
pub fn parse_exons<R: BufRead>(
    reader: &mut GtfReader<R>,
    genome_id_map: &mut HashMap<Vec<u8>, usize>,
    gene_id_map: &mut HashMap<Vec<u8>, usize>,
    transcript_id_map: &mut HashMap<Vec<u8>, usize>,
) -> Result<HashMap<usize, Vec<ExonRecord>>> {
    let mut transcript_records = HashMap::new();
    while let Some(Ok(record)) = reader.next() {
        if record.feature == b"exon" {
            let exon_record =
                ExonRecord::from_gtf_record(record, genome_id_map, gene_id_map, transcript_id_map)?;
            transcript_records
                .entry(exon_record.transcript())
                .or_insert_with(Vec::new)
                .push(exon_record);
        }
    }
    Ok(transcript_records)
}
