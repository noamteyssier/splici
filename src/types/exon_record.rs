use anyhow::{bail, Result};
use bedrs::Strand;
use gtftools::GtfRecord;
use hashbrown::HashMap;

use crate::Giv;

use super::Translater;

#[derive(Debug, Clone, Copy, Hash)]
pub struct ExonRecord {
    genome_id: usize,
    start: usize,
    end: usize,
    gene_id: usize,
    transcript_id: usize,
    strand: Strand,
}
impl From<ExonRecord> for Giv {
    fn from(record: ExonRecord) -> Self {
        Giv::new(record.genome_id, record.start, record.end, record.strand)
    }
}
impl From<&ExonRecord> for Giv {
    fn from(record: &ExonRecord) -> Self {
        Giv::new(record.genome_id, record.start, record.end, record.strand)
    }
}
impl ExonRecord {
    pub fn new(
        genome_id: usize,
        start: usize,
        end: usize,
        gene_id: usize,
        transcript_id: usize,
        strand: Strand,
    ) -> Self {
        Self {
            genome_id,
            start,
            end,
            gene_id,
            transcript_id,
            strand,
        }
    }

    pub fn from_gtf_record(
        record: GtfRecord,
        genome_translater: &mut Translater,
        gene_translater: &mut Translater,
        transcript_translater: &mut Translater,
        transcript_to_gene: &mut HashMap<usize, usize>,
    ) -> Result<Self> {
        // insert gene_id into gene_map
        let gene_id = if let Some(gene_id) = record.attribute.gene_id {
            gene_id
        } else {
            bail!("Missing Gene ID");
        };

        let transcript_id = if let Some(transcript_id) = record.attribute.transcript_id {
            transcript_id
        } else {
            bail!("Missing Transcript ID");
        };

        if !gene_translater.has_name(&gene_id) {
            let gene_idx = gene_translater.len();
            gene_translater.insert(gene_id.clone(), gene_idx);
        }

        if !genome_translater.has_name(&record.seqname) {
            let genome_idx = genome_translater.len();
            genome_translater.insert(record.seqname.clone(), genome_idx);
        }

        if !transcript_translater.has_name(&transcript_id) {
            let transcript_idx = transcript_translater.len();
            transcript_translater.insert(transcript_id.clone(), transcript_idx);
        }

        // get genome_id and gene_id
        let genome_id = genome_translater.get_id(&record.seqname).unwrap();
        let gene_id = gene_translater.get_id(&gene_id).unwrap();
        let transcript_id = transcript_translater.get_id(&transcript_id).unwrap();

        // insert transcript_id and gene_id into transcript_to_gene
        transcript_to_gene.insert(*transcript_id, *gene_id);

        let start = record.start;
        let end = record.end;
        let strand = match record.strand[0] {
            b'+' => Strand::Forward,
            b'-' => Strand::Reverse,
            _ => bail!("Invalid strand"),
        };

        // build GeneRecord
        Ok(Self::new(
            *genome_id,
            start,
            end,
            *gene_id,
            *transcript_id,
            strand,
        ))
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
