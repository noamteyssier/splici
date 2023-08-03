use crate::{types::ExonRecord, Giv, GivSet, IdMap, NameMap};
use anyhow::{bail, Result};
use bedrs::{Container, Coordinates, Internal, Merge};
use gtftools::GtfReader;
use hashbrown::HashMap;
use noodles::core::{Position, Region};
use std::{hash::Hash, io::BufRead, str::from_utf8};

/// Flips the K, V pairs in a `HashMap`
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
pub fn interval_to_region(giv: &Giv, genome_name_map: &IdMap) -> Result<Region> {
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

pub fn build_exon_set(exon_set: &[ExonRecord]) -> GivSet {
    let mut exon_set = exon_set
        .iter()
        .map(std::convert::Into::into)
        .collect::<GivSet>();
    exon_set.sort();
    exon_set
}

pub fn build_interval_set(vec: &[Giv]) -> GivSet {
    let mut set = GivSet::from_iter(vec.iter().copied());
    set.sort();
    set
}

pub fn get_introns(giv_set: GivSet, extension: Option<usize>) -> Result<Vec<Giv>> {
    let introns = giv_set
        .internal()?
        .map(|mut intron| {
            if let Some(extension) = extension {
                intron.extend(&extension);
            }
            intron
        })
        .collect::<Vec<Giv>>();
    Ok(introns)
}

pub fn merge_interval_set(giv_set: GivSet) -> GivSet {
    giv_set
        .merge()
        .unwrap()
        .intervals()
        .iter()
        .copied()
        .collect()
}

/// Parse the exons from the GTF file and return a `HashMap` of transcript id to
/// a Vec of `ExonRecords`
pub fn parse_exons<R: BufRead>(
    reader: &mut GtfReader<R>,
    genome_id_map: &mut NameMap,
    gene_id_map: &mut NameMap,
    transcript_id_map: &mut NameMap,
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
