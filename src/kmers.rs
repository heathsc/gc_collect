use serde::Deserialize;

use crate::{
    cli::Config,
    kmcv::{Kmcv, KmcvHeaderCore},
};

pub type KmerType = u32;

#[derive(Clone, Deserialize)]
pub struct KmerCounts {
    kmcv: KmcvHeaderCore,
    total_reads: u32,
    mapped_reads: u32,
    total_bases: u64,
    mapped_bases: u64,
    counts: Vec<(u32, u64)>,
}

impl KmerCounts {
    pub fn kmer_coverage(&self, cfg: &Config) -> anyhow::Result<Option<KmerCoverage>> {
        if let Some(kmcv) = cfg.kmcv() {
            self.get_coverage_stats(kmcv)
        } else {
            warn!("Cannot process kmer coverage without an input kmer file (use -k option)");
            Ok(None)
        }
    }

    fn get_coverage_stats(&self, kmcv: &Kmcv) -> anyhow::Result<Option<KmerCoverage>> {
        Ok(None)
    }

    pub fn add(&mut self, other: &Self) -> anyhow::Result<()> {
        if self.kmcv != other.kmcv {
            Err(anyhow!(
                "Cannot merge datasets as kmer files are not compatible"
            ))
        } else {
            self.total_reads += other.total_reads;
            self.total_bases += other.total_bases;
            self.mapped_reads += other.mapped_reads;
            self.mapped_bases += other.mapped_bases;

            assert_eq!(self.counts.len(), other.counts.len());
            for (p, q) in self.counts.iter_mut().zip(other.counts.iter()) {
                p.1 += q.1
            }

            Ok(())
        }
    }
}

#[derive(Debug)]
pub struct KmerCoverage {}
