use std::io::BufRead;

use crate::kmers::KmerType;
use anyhow::Context;
use log::{log_enabled, Level::Trace};
use serde::Deserialize;

fn get_u16_from_slice(p: &[u8]) -> u16 {
    u16::from_le_bytes(p.try_into().expect("Slice has wrong size"))
}

fn get_u32_from_slice(p: &[u8]) -> u32 {
    u32::from_le_bytes(p.try_into().expect("Slice has wrong size"))
}

#[derive(Clone, Debug, Eq, PartialEq, Deserialize)]
pub struct KmcvHeaderCore {
    version: [u8; 2],
    kmer_length: u8,
    max_hits: u8,
    n_contigs: u32,
    n_targets: u32,
    rnd_id: u32,
}

#[derive(Clone, Debug)]
pub struct KmcvHeader {
    core: KmcvHeaderCore,
}

impl KmcvHeader {
    pub fn read<R: BufRead>(rdr: &mut R) -> anyhow::Result<Self> {
        let mut buf = [0; 52];
        rdr.read_exact(&mut buf)
            .with_context(|| "Error reading header from kmer file")?;

        if buf[0..4] != [b'K', b'M', b'C', b'V'] {
            return Err(anyhow!(
                "Incorrect magic number from header block of  kmer file"
            ));
        }
        let version = [buf[4], buf[5]];
        if version[0] != 2 {
            return Err(anyhow!("Incorrect version for kmer file (expected V2)"));
        }
        let kmer_length = buf[6];
        let max_hits = buf[7];
        if (kmer_length as u32) << 1 > KmerType::BITS {
            return Err(anyhow!("Kmer length {kmer_length} too large for KmerType"));
        }
        let rnd_id = get_u32_from_slice(&buf[8..12]);
        let n_contigs = get_u32_from_slice(&buf[12..16]);
        let n_targets = get_u32_from_slice(&buf[16..20]);

        Ok(Self {
            core: KmcvHeaderCore {
                version,
                kmer_length,
                max_hits,
                rnd_id,
                n_contigs,
                n_targets,
            },
        })
    }
}

pub struct Target {
    start: u32,
    end: u32,
}

impl Target {
    #[inline]
    pub fn size(&self) -> u32 {
        self.end + 1 - self.start
    }
}
impl Target {
    fn read<R: BufRead>(rdr: &mut R, n_contigs: u32) -> anyhow::Result<(Self, u32)> {
        let mut buf = [0u8; 12];
        rdr.read_exact(&mut buf)
            .with_context(|| "Error reading target block from kmer file")?;

        let contig = Self::get_contig(&buf[..4], n_contigs)?;

        let (start, end) = Self::get_start_end(&buf[4..])?;

        Ok((Self { start, end }, contig))
    }

    fn get_contig(buf: &[u8], n_contigs: u32) -> anyhow::Result<u32> {
        let contig = get_u32_from_slice(&buf[..4]);
        if contig >= n_contigs {
            Err(anyhow!(
                "Contig id {contig}  from target definition not in range"
            ))
        } else {
            Ok(contig)
        }
    }

    fn get_start_end(buf: &[u8]) -> anyhow::Result<(u32, u32)> {
        let start = get_u32_from_slice(&buf[..4]);
        let end = get_u32_from_slice(&buf[4..]);
        if end < start {
            Err(anyhow!("End coordinate of target less than start"))
        } else {
            Ok((start, end))
        }
    }
}

pub struct KContig {
    name: Box<str>,
    targets: Vec<u32>,
}

impl KContig {
    fn read<R: BufRead>(rdr: &mut R) -> anyhow::Result<Self> {
        let mut buf = [0u8; 2];

        rdr.read_exact(&mut buf)
            .with_context(|| "Error reading string length from kmer file")?;

        let l = get_u16_from_slice(&buf) as usize;
        if l == 0 {
            return Err(anyhow!("Contig name length is zero"));
        }

        let mut s = String::with_capacity(l);
        while s.len() < l {
            let p = rdr
                .fill_buf()
                .with_context(|| "Error while reading contig name")?;
            let m = l - s.len();
            let n = m.min(p.len());
            let s1 = std::str::from_utf8(&p[..n]).with_context(|| "Contig name not utf8")?;
            s.push_str(s1);
            rdr.consume(n);
            if n == m {
                break;
            }
        }
        trace!("Read contig {s}");
        let name = s.into_boxed_str();
        let targets = Vec::new();
        Ok(Self { name, targets })
    }
}

pub struct Kmcv {
    header: KmcvHeader,
    contigs: Vec<KContig>,
    targets: Vec<Target>,
}

impl Kmcv {
    pub fn read<R: BufRead>(rdr: &mut R) -> anyhow::Result<Self> {
        debug!("Reading header from kmer file");
        let header =
            KmcvHeader::read(rdr).with_context(|| "Error reading header from kmer file")?;

        let n_ctgs = header.core.n_contigs as usize;
        let n_targets = header.core.n_targets as usize;

        let mut kmcv = Self {
            header,
            contigs: Vec::with_capacity(n_ctgs),
            targets: Vec::with_capacity(n_targets),
        };

        debug!("Reading contig blocks from kmer file");
        kmcv.read_contig_blocks(rdr)
            .with_context(|| "Error reading contig information")?;

        debug!("Reading target blocks from kmer file");
        kmcv.read_target_blocks(rdr)
            .with_context(|| "Error reading target information")?;

        Ok(kmcv)
    }

    pub fn get_target_size(&self, ix: usize) -> Option<u32> {
        self.targets.get(ix).map(|t| t.size())
    }

    /// Private functions
    fn read_contig_blocks<R: BufRead>(&mut self, rdr: &mut R) -> anyhow::Result<()> {
        self.contigs.clear();
        for _ in 0..self.header.core.n_contigs {
            self.contigs.push(KContig::read(rdr)?)
        }
        Ok(())
    }

    fn read_target_blocks<R: BufRead>(&mut self, rdr: &mut R) -> anyhow::Result<()> {
        let n_contigs = self.header.core.n_contigs;

        for ix in 0..self.header.core.n_targets {
            let (target, contig) = Target::read(rdr, n_contigs)?;
            self.contigs[contig as usize].targets.push(ix);
            self.targets.push(target);
        }

        if log_enabled!(Trace) {
            for ctg in self.contigs.iter() {
                debug!(
                    "Contig {} number of targets: {}",
                    ctg.name,
                    ctg.targets.len()
                )
            }
        }

        Ok(())
    }
}
