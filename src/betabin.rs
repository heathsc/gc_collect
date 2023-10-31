use anyhow::Context;
use std::{fs::File, io::Write, path::Path};

use compress_io::compress::CompressIo;
use libm::lgamma;

use crate::{
    gauss_legendre::gauss_legendre_64,
    reference::{GcHistKey, GcHistVal},
};

pub fn lbeta(a: f64, b: f64) -> f64 {
    lgamma(a) + lgamma(b) - lgamma(a + b)
}

pub fn mean_gc(cts: &[(GcHistKey, GcHistVal)]) -> f64 {
    let mut ct = [0.0; 2];
    for (a, b) in cts.iter().map(|(k, v)| {
        let (a, b) = k.counts();
        let w = v.count();
        (a * w, b * w)
    }) {
        ct[0] += a;
        ct[1] += b;
    }
    assert!(ct[0] + ct[1] > 0.0);
    ct[1] / (ct[0] + ct[1])
}

fn prob_func(x: f64, cts: &[(GcHistKey, GcHistVal)]) -> f64 {
    let lnx = x.ln();
    let lnx1 = (1.0 - x).ln();
    let (l, tot) = cts.iter().fold((0.0, 0.0), |(l, t), (c, v)| {
        let (a, b) = c.counts();
        let z = v.count();
        (l + (lnx * b + lnx1 * a - v.beta_a_b()).exp() * z, t + z)
    });
    l / tot
}
fn kl_distance_func(
    x: f64,
    cts: &[(GcHistKey, GcHistVal)],
    ref_dist: &[(GcHistKey, GcHistVal)],
) -> f64 {
    assert!(x > 0.0 && x < 1.0);
    let p = prob_func(x, cts);
    let q = prob_func(x, ref_dist);
    p * (p / q).ln()
}

pub fn kl_distance(cts: &[(GcHistKey, GcHistVal)], ref_dist: &[(GcHistKey, GcHistVal)]) -> f64 {
    gauss_legendre_64(|x| kl_distance_func(x, cts, ref_dist), 0.0, 1.0)
}

const GC_HIST_BINS: usize = 1000;

pub fn output_gc_hist(
    path: &Path,
    cts: &[(GcHistKey, GcHistVal)],
    ref_cts: Option<&[(GcHistKey, GcHistVal)]>,
) -> anyhow::Result<()> {
    let mut path1 = path.to_path_buf();
    path1.set_extension("gc_hist.tsv");

    let mut wrt = CompressIo::new()
        .path(&path1)
        .bufwriter()
        .with_context(|| "Could not open output gc distribution file")?;

    let mut lnp = Vec::with_capacity(GC_HIST_BINS);
    let mut tmp = Vec::with_capacity(GC_HIST_BINS);

    let bin_width = 1.0 / (GC_HIST_BINS as f64);

    for i in 0..GC_HIST_BINS {
        let x = bin_width * (0.5 + (i as f64));
        lnp.push((x, x.ln(), (1.0 - x).ln()))
    }

    let contrib = |c: &[(GcHistKey, GcHistVal)], tmp: &mut Vec<f64>, h: &mut [f64]| {
        tmp.clear();
        let mut t = 0.0;
        for (b, a, v) in c.iter().map(|(key, v)| {
            let (r, s) = key.counts();
            (r, s, v)
        }) {
            let x = v.count();
            t += x;
            let konst = v.beta_a_b();
            tmp.clear();
            let mut z = 0.0;
            for (_, lnp, lnp1) in lnp.iter() {
                let p = (lnp * a + lnp1 * b - konst).exp();
                z += p;
                tmp.push(p);
            }
            for (p, q) in tmp.iter().zip(h.iter_mut()) {
                *q += x * p / z
            }
        }
        t
    };

    let mut hist = vec![0.0; GC_HIST_BINS];
    let t = contrib(cts, &mut tmp, &mut hist);

    let rhist = ref_cts.map(|r| {
        let mut h = vec![0.0; GC_HIST_BINS];
        let t = contrib(r, &mut tmp, &mut h);
        (h, t)
    });

    write!(wrt, "GC\tSample")?;
    if rhist.is_some() {
        write!(wrt, "\tReference")?
    }
    writeln!(wrt)?;
    let z = GC_HIST_BINS as f64;
    for i in 0..1000 {
        write!(wrt, "{}\t{}", lnp[i].0, hist[i] * z / t)?;
        if let Some((rh, t1)) = rhist.as_ref() {
            write!(wrt, "\t{}", rh[i] * z / t1)?;
        }
        writeln!(wrt)?
    }

    Ok(())
}
