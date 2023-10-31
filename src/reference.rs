use std::{collections::HashMap, fs::File, io::BufReader, path::Path};

use serde::Deserialize;

use anyhow::Context;
use compress_io::compress::CompressIo;
use serde_json::{from_reader, from_value, Value};

use crate::betabin::lbeta;

fn get_value<'a>(js: &'a Value, ix: &str) -> anyhow::Result<&'a Value> {
    js.get(ix)
        .ok_or_else(|| anyhow!("Could not find {ix} field"))
}

#[derive(Deserialize)]
struct RSCounts {
    counts: HashMap<String, u64>,
    bisulfite_counts: Option<HashMap<String, u64>>,
}

#[derive(Deserialize)]
struct RawRef {
    read_lengths: Vec<u32>,
    read_length_specific_counts: HashMap<u32, RSCounts>,
}
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct GcHistKey(f64, f64);

impl GcHistKey {
    pub fn counts(&self) -> (f64, f64) {
        (self.0, self.1)
    }
}

impl GcHistKey {
    pub fn from_str(s: &str) -> anyhow::Result<Self> {
        if let Some((s1, s2)) = s.split_once(':') {
            let c1 = s1.parse::<u32>()?;
            let c2 = s2.parse::<u32>()?;
            Ok(Self(c1 as f64, c2 as f64))
        } else {
            Err(anyhow!("counts keys not in correct format"))
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct GcHistVal {
    count: f64,
    beta_a_b: f64,
}

impl GcHistVal {
    pub fn make(k: &GcHistKey, c: u64) -> Self {
        let (a, b) = k.counts();
        Self {
            count: c as f64,
            beta_a_b: lbeta(a + 1.0, b + 1.0),
        }
    }

    pub fn count(&self) -> f64 {
        self.count
    }
    pub fn beta_a_b(&self) -> f64 {
        self.beta_a_b
    }
}

pub struct Counts {
    regular: Vec<(GcHistKey, GcHistVal)>,
    bisulfite: Option<Vec<(GcHistKey, GcHistVal)>>,
}

impl Counts {
    fn from_rs_counts(mut rs: RSCounts) -> anyhow::Result<Self> {
        let RSCounts {
            mut counts,
            mut bisulfite_counts,
        } = rs;

        let make = |k: &String, v| -> anyhow::Result<_> {
            let key = GcHistKey::from_str(k)?;
            let val = GcHistVal::make(&key, v);
            Ok((key, val))
        };

        let mut regular = Vec::with_capacity(counts.len());
        for (k, v) in counts.drain() {
            let (key, val) = make(&k, v)?;
            regular.push((key, val));
        }
        let bisulfite = match bisulfite_counts.take() {
            Some(mut h) => {
                let mut b = Vec::with_capacity(h.len());
                for (k, v) in h.drain() {
                    let (key, val) = make(&k, v)?;
                    b.push((key, val));
                }
                Some(b)
            }
            None => None,
        };
        Ok(Self { regular, bisulfite })
    }

    pub fn regular(&self) -> &[(GcHistKey, GcHistVal)] {
        &self.regular
    }
    pub fn bisulfite(&self) -> Option<&[(GcHistKey, GcHistVal)]> {
        self.bisulfite.as_deref()
    }
}
pub struct RefDist {
    read_lengths: Vec<u32>,
    read_length_specific_counts: HashMap<u32, Counts>,
}

impl RefDist {
    fn from_raw(mut raw: RawRef) -> anyhow::Result<Self> {
        let RawRef {
            read_lengths,
            read_length_specific_counts: mut rlsc,
        } = raw;
        let mut read_length_specific_counts = HashMap::with_capacity(rlsc.len());
        for (k, v) in rlsc.drain() {
            read_length_specific_counts.insert(k, Counts::from_rs_counts(v)?);
        }
        Ok(Self {
            read_lengths,
            read_length_specific_counts,
        })
    }

    pub fn from_json_file<P: AsRef<Path>>(p: P) -> anyhow::Result<Self> {
        let p = p.as_ref();
        let rdr = CompressIo::new()
            .path(p)
            .bufreader()
            .with_context(|| format!("Could not open {} for input", p.display()))?;

        let mut raw: RawRef =
            from_reader(rdr).with_context(|| format!("Error parsing JSON file {}", p.display()))?;

        info!("Reference distributions read from {}", p.display());

        Self::from_raw(raw)
    }

    pub fn get_closest_reference(&self, rl: u32) -> (u32, &Counts) {
        let rlens = &self.read_lengths;
        let closest_ix = rlens[1..].iter().enumerate().fold(0, |k, (i, l)| {
            if rl.abs_diff(*l) < rl.abs_diff(rlens[k]) {
                i + 1
            } else {
                k
            }
        });
        let rl1 = rlens[closest_ix];
        (rl1, &self.read_length_specific_counts[&rl1])
    }
}
