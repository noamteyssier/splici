use anyhow::{bail, Result};
use bedrs::GenomicInterval;
use gtftools::GtfRecord;
use hashbrown::HashMap;

#[derive(Debug, Clone, Copy, Hash)]
pub struct ExonRecord {
    genome_id: usize,
    start: usize,
    end: usize,
    gene_id: usize,
    transcript_id: usize,
}
impl Into<GenomicInterval<usize>> for ExonRecord {
    fn into(self) -> GenomicInterval<usize> {
        GenomicInterval::new(self.genome_id, self.start, self.end)
    }
}
impl Into<GenomicInterval<usize>> for &ExonRecord {
    fn into(self) -> GenomicInterval<usize> {
        GenomicInterval::new(self.genome_id, self.start, self.end)
    }
}
impl ExonRecord {
    pub fn new(
        genome_id: usize,
        start: usize,
        end: usize,
        gene_id: usize,
        transcript_id: usize,
    ) -> Self {
        Self {
            genome_id,
            start,
            end,
            gene_id,
            transcript_id,
        }
    }

    pub fn from_gtf_record(
        record: GtfRecord,
        genome_map: &mut HashMap<Vec<u8>, usize>,
        gene_map: &mut HashMap<Vec<u8>, usize>,
        transcript_map: &mut HashMap<Vec<u8>, usize>,
    ) -> Result<Self> {
        // insert gene_id into gene_map
        let gene_id = if let Some(gene_id) = record.attribute.gene_id {
            gene_id
        } else {
            bail!("Missing Gene ID");
        };
        if !gene_map.contains_key(&gene_id) {
            let gene_idx = gene_map.len();
            gene_map.insert(gene_id.clone(), gene_idx);
        }

        let transcript_id = if let Some(transcript_id) = record.attribute.transcript_id {
            transcript_id
        } else {
            bail!("Missing Transcript ID");
        };

        // insert genome_id into genome_map
        if !genome_map.contains_key(&record.seqname) {
            let genome_idx = genome_map.len();
            genome_map.insert(record.seqname.clone(), genome_idx);
        }

        // insert transcript_id into transcript_map
        if !transcript_map.contains_key(&transcript_id) {
            let transcript_idx = transcript_map.len();
            transcript_map.insert(transcript_id.clone(), transcript_idx);
        }

        // get genome_id and gene_id
        let genome_id = genome_map[&record.seqname];
        let gene_id = gene_map[&gene_id];
        let transcript_id = transcript_map[&transcript_id];
        let start = record.start;
        let end = record.end;

        // build GeneRecord
        Ok(Self::new(genome_id, start, end, gene_id, transcript_id))
    }

    pub fn genome(&self) -> usize {
        self.genome_id
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }

    pub fn gene(&self) -> usize {
        self.gene_id
    }

    pub fn transcript(&self) -> usize {
        self.transcript_id
    }
}
