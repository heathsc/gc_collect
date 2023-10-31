use anyhow::Context;
use compress_io::compress::CompressIo;
use crossbeam_channel::{Sender, Receiver};
use std::{io::Write, path::Path};

use crate::{
    betabin::*,
    cli::Config,
    read::{read_json, BisulfiteType, DataSet},
    simple_regression::*,
};

#[derive(Debug)]
pub struct DataResults {
    mean_gc: f64,
    ref_mean_gc: Option<f64>,
    kl_distance: Option<f64>,
    regression: Option<Vec<SimpleRegression>>,
}

fn compare_to_reference(
    cfg: &Config,
    path: &Path,
    d: &DataSet,
) -> anyhow::Result<(Option<f64>, Option<f64>)> {
    let (r, kl, gc) = match cfg.ref_dist() {
        Some(r) => {
            let (rl, counts) = r.get_closest_reference(d.max_read_len() as u32);
            trace!(
                "Using reference length {rl} for actual length {}",
                d.max_read_len()
            );

            let ref_counts = match d.bisulfite() {
                BisulfiteType::None => Some(counts.regular()),
                _ => counts.bisulfite(),
            };

            (
                ref_counts,
                ref_counts.map(|ref_counts| kl_distance(d.gc_counts(), ref_counts)),
                ref_counts.map(mean_gc),
            )
        }
        None => (None, None, None),
    };

    output_gc_hist(path, d.gc_counts(), r).with_context(|| "Error writing gc distribution file")?;
    Ok((kl, gc))
}

fn base_content_regressions(d: &DataSet) -> Option<Vec<SimpleRegression>> {
    let ct = d.per_pos_cts();
    let l = ct.len();
    let x0 = l / 3;
    if l - x0 < 3 {
        return None;
    }
    let scale = (l - x0) as f64;
    let mut obs = Vec::with_capacity(l - x0);
    let mut res = Vec::with_capacity(4);
    for ix in 0..4 {
        obs.clear();
        for (x, y) in ct[x0..]
            .iter()
            .map(|c| {
                let s = c.cts()[..4].iter().sum::<u64>() as f64;
                c.cts()[ix] as f64 / s
            })
            .enumerate()
        {
            obs.push(((x as f64) / scale, y))
        }
        let reg = match simple_regression(&obs) {
            Ok(r) => r,
            Err(e) => {
                warn!("Could not perform regression: {:?}", e);
                return None;
            }
        };

        res.push(reg)
    }
    Some(res)
}

fn output_per_cycle_bases(d: &DataSet, p: &Path) -> anyhow::Result<()> {
    let mut path = p.to_path_buf();
    path.set_extension("base_dist.tsv");
    let mut wrt = CompressIo::new()
        .path(&path)
        .bufwriter()
        .with_context(|| "Could not open output file")?;

    writeln!(wrt, "Cycle\tA\tC\tG\tT")?;
    let trim = d.trim();
    let cts = d.per_pos_cts();
    for (i, ct) in cts.iter().enumerate() {
        let s = ct.cts()[..4].iter().sum::<u64>() as f64;
        write!(wrt, "{}", i + 1 + trim)?;
        for k in [0, 1, 3, 2] {
            let y = (ct.cts()[k] as f64) / s;
            write!(wrt, "\t{:.5}", y)?;
        }
        writeln!(wrt)?
    }
    Ok(())
}

fn analyze_dataset(cfg: &Config, path: &Path, d: &DataSet) -> anyhow::Result<DataResults> {
    output_per_cycle_bases(d, path).with_context(|| "Error writing per cycle base distribution")?;
    let mean_gc = mean_gc(d.gc_counts());
    let (kl_distance, ref_mean_gc) = compare_to_reference(cfg, path, d)?;
    let regression = base_content_regressions(d);
    Ok(DataResults {
        mean_gc,
        kl_distance,
        ref_mean_gc,
        regression,
    })
}
fn process_file(cfg: &Config, p: &Path) -> anyhow::Result<DataResults> {
    trace!("Reading from {}", p.display());
    let d = read_json(p).with_context(|| format!("Error reading from {}", p.display()))?;
    let dres = analyze_dataset(cfg, p, &d)?;
    /* print!(
        "{:?}\t{:?}\t{:?}",
        dres.mean_gc, dres.ref_mean_gc, dres.kl_distance
    );
    if let Some(v) = dres.regression.as_ref() {
        for r in v {
            print!("\t{:?}\t{:?}", r.slope().estimate(), r.slope().p())
        }
    }
    println!(); */
    Ok(dres)
}

pub fn process_thread(cfg: &Config, ix: usize, rx: Receiver<&Path>, sd: Sender<DataResults>) -> anyhow::Result<()> {
    debug!("Process thread {ix} starting up");
    while let Ok(p) = rx.recv() {
        trace!(
            "Process thread {ix} received file {} for processing",
            p.display()
        );
        let dres = process_file(cfg, p)?;
        trace!("Process thread {ix} finished processing file {}", p.display());
        sd.send(dres).with_context(|| "Error sending results to output thread")?
    }
    debug!("Process thread {ix} closing down");
    Ok(())
}
