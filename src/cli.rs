use std::path::{Path, PathBuf};

use anyhow::Context;
use compress_io::compress::CompressIo;

mod cli_model;

use crate::{kmcv::Kmcv, reference::RefDist};
pub use cli_model::MergeKey;

pub struct Config {
    input_files: Vec<PathBuf>,
    output_file: Option<PathBuf>,
    ref_dist: Option<RefDist>,
    threads: usize,
    regression: bool,
    kmcv: Option<Kmcv>,
    merge_key: Option<MergeKey>,
}

impl Config {
    pub fn input_files(&self) -> &[PathBuf] {
        &self.input_files
    }
    pub fn output_file(&self) -> Option<&Path> {
        self.output_file.as_deref()
    }
    pub fn threads(&self) -> usize {
        self.threads
    }
    pub fn ref_dist(&self) -> Option<&RefDist> {
        self.ref_dist.as_ref()
    }
    pub fn regression(&self) -> bool {
        self.regression
    }
    pub fn kmcv(&self) -> Option<&Kmcv> {
        self.kmcv.as_ref()
    }
    pub fn merge_key(&self) -> Option<MergeKey> {
        self.merge_key
    }
}
pub fn handle_cli() -> anyhow::Result<Config> {
    let c = cli_model::cli_model();
    let m = c.get_matches();
    super::utils::init_log(&m);

    let input_files: Vec<PathBuf> = m
        .get_many("input")
        .expect("Missing required input argument")
        .map(|p: &PathBuf| p.to_owned())
        .collect();

    let output_file = m.get_one::<PathBuf>("output").map(|p| p.to_owned());
    let threads = m
        .get_one::<u64>("threads")
        .map(|x| *x as usize)
        .unwrap_or_else(|| num_cpus::get().min(input_files.len()));

    let ref_dist = match m.get_one::<PathBuf>("ref") {
        Some(p) => Some(RefDist::from_json_file(&p).with_context(|| {
            format!(
                "Error reading reference distributions from JSON file {}",
                p.display()
            )
        })?),
        None => None,
    };

    let regression = m.get_flag("regression");
    
    let merge_key = m.get_one::<MergeKey>("merge_by").copied().or_else(|| {
        if m.get_flag("merge") {
            Some(MergeKey::Default)
        } else {
            None
        }
    });

    let kmcv = match m.get_one::<PathBuf>("kmers") {
        Some(p) => {
            let mut rdr = CompressIo::new()
                .path(p)
                .bufreader()
                .with_context(|| "Could not open kmer file for input")?;

            debug!("Opened kmer file for input");
            Some(
                Kmcv::read(&mut rdr)
                    .with_context(|| format!("Could not read kmer file {}", p.display()))?,
            )
        }
        None => None,
    };

    Ok(Config {
        input_files,
        output_file,
        merge_key,
        threads,
        ref_dist,
        regression,
        kmcv,
    })
}
