use anyhow::Context;
use std::path::{Path, PathBuf};

mod cli_model;

use crate::reference::RefDist;

pub struct Config {
    input_files: Vec<PathBuf>,
    output_file: Option<PathBuf>,
    ref_dist: Option<RefDist>,
    threads: usize,
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

    Ok(Config {
        input_files,
        output_file,
        threads,
        ref_dist,
    })
}
