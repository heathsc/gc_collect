use anyhow::Context;
use crossbeam_channel::{Receiver, Sender};
use std::{
    collections::{hash_map, HashMap},
    path::{Path, PathBuf},
};

use crate::{
    cli::{Config, MergeKey},
    read::{read_json, DataSet, Fli},
};

fn get_merge_key(fli: &mut Fli, mut m: MergeKey) -> anyhow::Result<(MergeKey, String)> {
    if matches!(m, MergeKey::Default) {
        m = fli
            .find_merge_key()
            .ok_or(anyhow!("Couldn't determine merge key type for dataset"))?
    }
    let key = fli
        .get_key(m)
        .ok_or(anyhow!("Couldn't establish merge key not dataset"))?;

    Ok((m, key))
}

fn merge_dataset(
    mut d: DataSet,
    m: MergeKey,
    hash: &mut HashMap<String, DataSet>,
) -> anyhow::Result<MergeKey> {
    let (m, key) = get_merge_key(d.fli_mut(), m)?;

    let path = PathBuf::from(&key);
    match hash.entry(key) {
        hash_map::Entry::Occupied(mut e) => e.get_mut().merge(&d)?,
        hash_map::Entry::Vacant(e) => {
            d.set_path(path);
            e.insert(d);
        }
    }

    Ok(m)
}

pub fn merge_thread(cfg: &Config, rx: Receiver<&Path>, sd: Sender<DataSet>) -> anyhow::Result<()> {
    debug!("Merge thread starting up");

    let mut merge_key = cfg.merge_key().expect("Cannot merge without a key!");

    let mut hash: HashMap<String, DataSet> = HashMap::new();

    while let Ok(p) = rx.recv() {
        trace!("Merge thread received file {} for reading", p.display());

        let d = read_json(p).with_context(|| format!("Error reading from {}", p.display()))?;
        merge_key = merge_dataset(d, merge_key, &mut hash)?;
    }

    debug!("Merge thread finished merging all input files. Sending results to process thread");

    for (_, mut d) in hash.drain() {
        d.mk_gc_counts()?;
        sd.send(d)
            .with_context(|| "Error sending results to process thread")?
    }

    debug!("Merge thread closing down");

    Ok(())
}
