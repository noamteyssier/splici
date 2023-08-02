use std::{fs::File, io::{BufReader, BufRead}, str::from_utf8};
use anyhow::Result;
use bedrs::{GenomicInterval, GenomicIntervalSet, Container, Merge, Internal, Coordinates};
use gtftools::GtfReader;
use hashbrown::HashMap;
use noodles::fasta::{fai, IndexedReader};
use crate::{utils::{flip_map, interval_to_region}, types::ExonRecord};

/// Parse the exons from the GTF file and return a HashMap of transcript id to
/// a Vec of `ExonRecords`
fn parse_exons<R: BufRead>(
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

fn populate_introns(
    transcript_records: &HashMap<usize, Vec<ExonRecord>>,
    extension: Option<usize>,
) -> HashMap<usize, Vec<GenomicInterval<usize>>> {
    transcript_records
        .values()
        .map(|exon_set| {
            // identify the parent gene
            let parent_id = exon_set[0].gene();

            // select the exon set of the transcript
            let mut exon_set = exon_set
                .iter()
                .map(|exon| exon.into())
                .collect::<GenomicIntervalSet<usize>>();
            exon_set.sort();

            // calculate the introns
            let intron_set = exon_set
                .internal()
                .expect("Could not calculate introns")
                .map(|mut intron| {
                    if let Some(ref ext) = extension {
                        intron.extend(ext);
                    }
                    intron
                });

            // return the parent gene id and the intron set
            (parent_id, intron_set)
        })
        // fold the intron sets into a HashMap of gene id to intron set
        .fold(
            HashMap::new(),
            |mut gene_introns, (parent_id, intron_set)| {
                intron_set.for_each(|intron| {
                    gene_introns
                        .entry(parent_id)
                        .or_insert_with(Vec::new)
                        .push(intron);
                });
                gene_introns
            },
        )
}

fn merge_introns(
    gene_introns: &HashMap<usize, Vec<GenomicInterval<usize>>>,
) -> HashMap<usize, GenomicIntervalSet<usize>> {
    gene_introns
        .into_iter()
        .map(|(gene_id, intron_vec)| {
            let intron_set = GenomicIntervalSet::from_unsorted(intron_vec.clone());
            let merged_introns = intron_set
                .merge()
                .expect("Could not merge introns")
                .intervals()
                .into_iter()
                .cloned()
                .collect::<GenomicIntervalSet<usize>>();
            (*gene_id, merged_introns)
        })
        .collect::<HashMap<usize, GenomicIntervalSet<usize>>>()
}

pub fn run(
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

    let mut genome_id_map = HashMap::new();
    let mut gene_id_map = HashMap::new();
    let mut transcript_id_map = HashMap::new();

    let transcript_records = parse_exons(
        &mut reader,
        &mut genome_id_map,
        &mut gene_id_map,
        &mut transcript_id_map,
    )?;

    let genome_name_map = flip_map(&genome_id_map);
    let gene_name_map = flip_map(&gene_id_map);
    let transcript_id_map = flip_map(&transcript_id_map);

    // build intron locations by subtracting exons from parent
    let gene_introns = populate_introns(&transcript_records, extension);
    let merged_introns = merge_introns(&gene_introns);

    // extract exon sequences
    for transcript in transcript_records.keys() {
        let full_transcript_name = transcript_id_map
            .get(transcript)
            .map(|x| from_utf8(x))
            .unwrap()?;
        let exon_set = transcript_records.get(transcript).unwrap();
        let mut exon_set = exon_set
            .iter()
            .map(|x| x.into())
            .collect::<GenomicIntervalSet<usize>>();
        exon_set.sort();
        let transcript_seq = exon_set
            .records()
            .iter()
            .map(|x| {
                let region = interval_to_region(x, &genome_name_map).unwrap();
                let query = fasta.query(&region).unwrap();
                let seq = query.sequence().as_ref();
                let seq = from_utf8(seq).unwrap();
                seq.to_string()

                // println!("{}\t{}\t{}\t{full_transcript_name}", region.name(), x.start(), x.end());
                // println!("{full_transcript_name} {start} {end} {seq}", start=x.start(), end=x.end());

                // println!(">{}", full_transcript_name);
                // println!("{}", seq);
            })
            .fold(String::new(), |mut acc, x| {
                acc.push_str(&x);
                acc
            });
        print!(">{full_transcript_name}\n{transcript_seq}\n");
    }

    // extract intron sequences
    for gene in merged_introns.keys() {
        let intron_set = merged_introns.get(gene).unwrap();
        let full_gene_name = gene_name_map.get(gene).unwrap();
        let full_gene_name = from_utf8(full_gene_name)?;
        intron_set
            .records()
            .iter()
            .enumerate()
            .for_each(|(idx, x)| {
                let region = interval_to_region(x, &genome_name_map).unwrap();
                let query = fasta.query(&region).unwrap();
                let seq = from_utf8(query.sequence().as_ref()).unwrap();
                print!(">{full_gene_name}-I.{idx}\n{seq}\n");
            })
    }
    Ok(())
}
