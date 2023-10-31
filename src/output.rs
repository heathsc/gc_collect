use std::io::Write;
use crossbeam_channel::Receiver;
use anyhow::Result;
use compress_io::compress::CompressIo;

use crate::{
    cli::Config,
    process::DataResults,
};

fn output_result(res: DataResults) -> anyhow::Result<()> {

    Ok(())
}

pub fn output_thread(cfg: &Config, rx: Receiver<DataResults>) -> anyhow::Result<()> {
    debug!("Output thread starting up");

    while let Ok(res) = rx.recv() {
        output_result(res)?
    }
    
    debug!("Output thread closing down");
    Ok(())
}