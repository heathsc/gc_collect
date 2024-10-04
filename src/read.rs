use std::{
    collections::{BTreeMap, HashMap},
    ffi::OsStr,
    fmt,
    path::{Path, PathBuf},
};

use anyhow::Context;
use compress_io::compress::CompressIo;
use serde::Deserialize;
use serde_json::from_reader;

use crate::{
    cli::MergeKey,
    kmers::KmerCounts,
    reference::{GcHistKey, GcHistVal},
};

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

#[derive(Clone, Debug, Deserialize)]
pub struct Fli {
    sample: Option<String>,
    barcode: Option<String>,
    library: Option<String>,
    flowcell: Option<String>,
    index: Option<String>,
    lane: Option<u8>,
    read_end: Option<u8>,
}

impl Fli {
    pub fn get_key(&self, key: MergeKey) -> Option<String> {
        match key {
            MergeKey::Sample => self.sample.as_ref().map(|x| x.to_owned()),
            MergeKey::Barcode => self.barcode.as_ref().map(|x| x.to_owned()),
            MergeKey::Library => self.library.as_ref().map(|x| x.to_owned()),
            MergeKey::Fli => self.fli(),
            MergeKey::Default => None,
        }
    }

    pub fn find_merge_key(&self) -> Option<MergeKey> {
        if self.sample.is_some() {
            Some(MergeKey::Sample)
        } else if self.barcode.is_some() {
            Some(MergeKey::Barcode)
        } else if self.library.is_some() {
            Some(MergeKey::Library)
        } else if self.flowcell.is_some() && self.lane.is_some() && self.index.is_some() {
            Some(MergeKey::Fli)
        } else {
            None
        }
    }

    pub fn fli(&self) -> Option<String> {
        if let (Some(fc), Some(lane), Some(index)) =
            (self.flowcell.as_ref(), self.lane, self.index.as_ref())
        {
            Some(format!("{}_{}_{}", fc, lane, index))
        } else {
            None
        }
    }

    fn find_common(&mut self, other: &Self) {
        if self.sample != other.sample {
            self.sample = None
        }
        if self.barcode != other.barcode {
            self.barcode = None
        }
        if self.library != other.library {
            self.library = None
        }
        if self.index != other.index {
            self.index = None
        }
        if self.lane != other.lane {
            self.lane = None
        }
        if self.flowcell != other.flowcell {
            self.flowcell = None
        }
        if self.read_end != other.read_end {
            self.read_end = None
        }
    }
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
        write!(f, "\t{}", self.barcode.as_deref().unwrap_or("NA"))?;
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

#[derive(Copy, Clone, Debug, Default)]
pub struct Counts([u64; 5]);

impl Counts {
    fn from_temp_counts(t: &TempCounts) -> Self {
        let n = t.N.unwrap_or(0) + t.Other.unwrap_or(0);
        let ct = Self([t.A, t.C, t.T, t.G, n]);
        ct
    }

    pub fn cts(&self) -> &[u64; 5] {
        &self.0
    }

    pub fn add(&mut self, other: &Self) {
        for i in 0..5 {
            self.0[i] += other.0[i];
        }
    }
}

#[derive(Deserialize)]
struct TempDataSet {
    trim: usize,
    min_qual: u8,
    max_read_length: usize,
    bisulfite: BisulfiteType,
    fli: Fli,
    cts: TempCounts,
    per_pos_cts: BTreeMap<u32, TempCounts>,
    gc_hash: HashMap<String, u64>,
    kmer_counts: Option<KmerCounts>,
}

#[derive(Clone)]
pub struct DataSet {
    path: PathBuf,
    trim: usize,
    min_qual: u8,
    max_read_length: usize,
    bisulfite: BisulfiteType,
    fli: Fli,
    cts: Counts,
    per_pos_cts: Vec<Counts>,
    gc_hash: HashMap<String, u64>,
    gc_counts: Option<Vec<(GcHistKey, GcHistVal)>>,
    kmer_counts: Option<KmerCounts>,
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
    pub fn gc_counts(&self) -> Option<&[(GcHistKey, GcHistVal)]> {
        self.gc_counts.as_deref()
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

    pub fn kmer_counts(&self) -> Option<&KmerCounts> {
        self.kmer_counts.as_ref()
    }

    pub fn fli_mut(&mut self) -> &mut Fli {
        &mut self.fli
    }

    pub fn path(&self) -> &Path {
        self.path.as_path()
    }

    pub fn set_path(&mut self, path: PathBuf) {
        self.path = path
    }

    pub fn mk_gc_counts(&mut self) -> anyhow::Result<()> {
        let mut gc_counts = Vec::with_capacity(self.gc_hash.len());
        for (k, v) in self.gc_hash.iter() {
            let key = GcHistKey::from_str(k)?;
            let val = GcHistVal::make(&key, *v);
            gc_counts.push((key, val));
        }
        self.gc_counts = Some(gc_counts);
        Ok(())
    }

    fn add_gc_hash(&mut self, other: &Self) {
        let gc_hash = &mut self.gc_hash;

        for (k, v) in other.gc_hash.iter() {
            match gc_hash.get_mut(k) {
                Some(x) => *x += *v,
                None => {
                    gc_hash.insert(k.clone(), *v);
                }
            }
        }
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
            gc_hash,
            kmer_counts,
        } = t;

        let cts = Counts::from_temp_counts(&tmp_cts);
        let l = tmp_ppc.len();
        assert!(max_read_length >= trim && max_read_length - trim == l);
        let mut per_pos_cts = Vec::with_capacity(l);
        for (ix, (k, v)) in tmp_ppc.iter().enumerate() {
            assert_eq!(*k as usize, ix + 1 + trim);
            per_pos_cts.push(Counts::from_temp_counts(v))
        }

        let s = OsStr::new("gz");
        let path = if p.extension() == Some(s) {
            PathBuf::from(p.file_stem().unwrap())
        } else {
            p.to_owned()
        };

        Ok(Self {
            path,
            trim,
            min_qual,
            max_read_length,
            bisulfite,
            fli,
            cts,
            per_pos_cts,
            gc_hash,
            gc_counts: None,
            kmer_counts,
        })
    }

    fn check_constants(&self, other: &Self) -> bool {
        self.trim == other.trim
            && self.min_qual == other.min_qual
            && self.bisulfite == other.bisulfite
            && ((self.kmer_counts.is_some() && other.kmer_counts.is_some())
                || (self.kmer_counts.is_none() && other.kmer_counts.is_none()))
    }

    fn add_counts(&mut self, other: &Self) {
        self.cts.add(&other.cts);
        self.per_pos_cts
            .resize_with(self.max_read_length, Default::default);
        for (c1, c2) in self.per_pos_cts.iter_mut().zip(other.per_pos_cts().iter()) {
            c1.add(c2)
        }
    }
    pub fn merge(&mut self, other: &Self) -> anyhow::Result<()> {
        if !self.check_constants(other) {
            Err(anyhow!(
                "Cannot merge datasets generated with different parameters"
            ))
        } else {
            self.max_read_length = self.max_read_length.max(other.max_read_length);
            self.fli.find_common(&other.fli);
            self.add_counts(other);
            self.add_gc_hash(other);
            if let Some(kc) = self.kmer_counts.as_mut() {
                kc.add(other.kmer_counts().as_ref().unwrap())?
            }
            Ok(())
        }
    }
}

pub fn read_json<P: AsRef<Path>>(p: P) -> anyhow::Result<DataSet> {
    let p = p.as_ref();

    let rdr = CompressIo::new()
        .path(p)
        .bufreader()
        .with_context(|| format!("Could not open {} for input", p.display()))?;
    let tmp: TempDataSet = from_reader(rdr).with_context(|| "Error parsing JSON file")?;
    DataSet::from_temp_dataset(tmp, p)
}
