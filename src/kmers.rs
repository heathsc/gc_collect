use std::fmt;

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
    pub fn kmer_coverage(&self, cfg: &Config) -> Option<KmerCoverage> {
        if let Some(kmcv) = cfg.kmcv() {
            Some(self.get_coverage(kmcv))
        } else {
            warn!("Cannot process kmer coverage without an input kmer file (use -k option)");
            None
        }
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

    fn get_coverage(&self, kmcv: &Kmcv) -> KmerCoverage {
        let mut v: Vec<_> = self
            .counts
            .iter()
            .enumerate()
            .map(|(target_ix, (_, bases))| {
                let target_size = kmcv.get_target_size(target_ix).expect("Bad target ix") as f64;
                *bases as f64 / target_size
            })
            .collect();

        v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let l = v.len();
        let mean = v.iter().sum::<f64>() / (l as f64);
        let quartiles = [v[l >> 2], v[l >> 1], v[(3 * l) >> 2]];
        KmerCoverage{mean, quartiles}    
    }
}

#[derive(Debug)]
pub struct KmerCoverage {
    mean: f64,
    quartiles: [f64; 3],
}

impl KmerCoverage {
    pub fn median(&self) -> f64 {
        self.quartiles[1]
    }
    
    pub fn iqr(&self) -> f64 {
        self.quartiles[2] - self.quartiles[0]
    }
    
    pub fn dispersion(&self) -> f64 {
        self.iqr() / (self.quartiles[0] + self.quartiles[1])
    }
}
impl fmt::Display for KmerCoverage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{:.6}\t{:.6}", self.mean, self.median(), self.median() / self.mean, self.dispersion())
    }
}
