use std::{
    collections::{BTreeMap, HashMap},
    fmt,
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use anyhow::Context;
use serde::Deserialize;
use serde_json::from_reader;

use crate::reference::{GcHistKey, GcHistVal};

#[derive(Debug, Clone, Copy, Eq, PartialEq, Deserialize)]
pub enum BisulfiteType {
    None = 0,
    Forward,
    Reverse,
    NonStranded,
}

impl fmt::Display for BisulfiteType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::None => "None",
                Self::Forward => "Forward",
                Self::Reverse => "Reverse",
                Self::NonStranded => "Non-stranded",
            }
        )
    }
}

#[derive(Debug, Deserialize)]
pub struct Fli {
    sample: Option<String>,
    library: Option<String>,
    flowcell: Option<String>,
    index: Option<String>,
    lane: Option<u8>,
    read_end: Option<u8>,
}

impl fmt::Display for Fli {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let output_opt_u8 = |x: Option<u8>, f: &mut fmt::Formatter| -> fmt::Result {
            if let Some(x) = x {
                write!(f, "\t{}", x)
            } else {
                write!(f, "\tNA")
            }
        };

        write!(f, "{}", self.sample.as_deref().unwrap_or("NA"))?;
        write!(f, "\t{}", self.library.as_deref().unwrap_or("NA"))?;
        write!(f, "\t{}", self.flowcell.as_deref().unwrap_or("NA"))?;
        write!(f, "\t{}", self.index.as_deref().unwrap_or("NA"))?;
        output_opt_u8(self.lane, f)?;
        output_opt_u8(self.read_end, f)
    }
}
#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
struct TempCounts {
    A: u64,
    C: u64,
    G: u64,
    T: u64,
    N: Option<u64>,
    Other: Option<u64>,
}

#[derive(Debug)]
pub struct Counts([u64; 5]);

impl Counts {
    fn from_temp_counts(t: &TempCounts) -> Self {
        let n = t.N.unwrap_or(0) + t.Other.unwrap_or(0);
        Self([t.A, t.C, t.T, t.G, n])
    }

    pub fn cts(&self) -> &[u64; 5] {
        &self.0
    }
}

#[derive(Debug, Deserialize)]
struct TempDataSet {
    trim: usize,
    min_qual: u8,
    max_read_length: usize,
    bisulfite: BisulfiteType,
    fli: Fli,
    cts: TempCounts,
    per_pos_cts: BTreeMap<u32, TempCounts>,
    gc_hash: HashMap<String, u64>,
}

pub struct DataSet {
    path: PathBuf,
    trim: usize,
    min_qual: u8,
    max_read_length: usize,
    bisulfite: BisulfiteType,
    fli: Fli,
    #[allow(unused)]
    cts: Counts,
    per_pos_cts: Vec<Counts>,
    gc_counts: Vec<(GcHistKey, GcHistVal)>,
}

impl fmt::Display for DataSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.fli,
            self.path.display(),
            self.bisulfite,
            self.trim,
            self.min_qual
        )?;
        Ok(())
    }
}

impl DataSet {
    pub fn gc_counts(&self) -> &[(GcHistKey, GcHistVal)] {
        &self.gc_counts
    }
    pub fn bisulfite(&self) -> &BisulfiteType {
        &self.bisulfite
    }
    pub fn max_read_len(&self) -> usize {
        self.max_read_length
    }
    pub fn trim(&self) -> usize {
        self.trim
    }

    pub fn per_pos_cts(&self) -> &[Counts] {
        &self.per_pos_cts
    }

    fn from_temp_dataset(t: TempDataSet, p: &Path) -> anyhow::Result<Self> {
        let TempDataSet {
            trim,
            min_qual,
            max_read_length,
            bisulfite,
            fli,
            cts: tmp_cts,
            per_pos_cts: tmp_ppc,
            gc_hash: mut tmp_hash,
        } = t;

        let cts = Counts::from_temp_counts(&tmp_cts);
        let l = tmp_ppc.len();
        assert!(max_read_length >= trim && max_read_length - trim == l);
        let mut per_pos_cts = Vec::with_capacity(l);
        for (ix, (k, v)) in tmp_ppc.iter().enumerate() {
            assert_eq!(*k as usize, ix + 1 + trim);
            per_pos_cts.push(Counts::from_temp_counts(v))
        }
        let mut gc_counts = Vec::with_capacity(tmp_hash.len());
        for (k, v) in tmp_hash.drain() {
            let key = GcHistKey::from_str(&k)?;
            let val = GcHistVal::make(&key, v);
            gc_counts.push((key, val));
        }
        let path = p.to_owned();

        Ok(Self {
            path,
            trim,
            min_qual,
            max_read_length,
            bisulfite,
            fli,
            cts,
            per_pos_cts,
            gc_counts,
        })
    }
}

pub fn read_json<P: AsRef<Path>>(p: P) -> anyhow::Result<DataSet> {
    let p = p.as_ref();
    let file =
        File::open(p).with_context(|| format!("Could not open {} for input", p.display()))?;

    let rdr = BufReader::new(file);
    let tmp: TempDataSet = from_reader(rdr).with_context(|| "Error parsing JSON file")?;
    DataSet::from_temp_dataset(tmp, p)
}
