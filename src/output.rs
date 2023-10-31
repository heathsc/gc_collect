use anyhow::Context;
use compress_io::compress::CompressIo;
use crossbeam_channel::Receiver;
use std::io::Write;

use crate::{cli::Config, process::DataResults, read::DataSet};

pub fn output_thread(cfg: &Config, rx: Receiver<(DataSet, DataResults)>) -> anyhow::Result<()> {
    debug!("Output thread starting up");

    let mut wrt = CompressIo::new()
        .opt_path(cfg.output_file())
        .bufwriter()
        .with_context(|| "Could not open output file")?;

    writeln!(
        wrt,
        "Sample\tLibrary\tFlowcell\tIndex\tLane\tRead-end\tFile\tBisulfite-type\tTrim\tMin-qual\tgc\tref-gc\tKL-distance\tb(A)\tlog10 p_b(A)\tb(C)\tlog10 p_b(C)\tb(G)\tlog10 p_b(G)\tb(T)\tlog10 p_b(T)")?;

    while let Ok((data, res)) = rx.recv() {
        writeln!(wrt, "{}\t{}", data, res)?
    }

    debug!("Output thread closing down");
    Ok(())
}
