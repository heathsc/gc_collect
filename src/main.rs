#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

use crossbeam_channel::{bounded, unbounded};
use crossbeam_utils::thread::{self, ScopedJoinHandle};

mod betabin;
mod cli;
mod gauss_legendre;
mod kmcv;
mod kmers;
mod merge;
mod output;
mod process;
mod read;
mod reference;
mod simple_regression;
mod utils;

use cli::Config;
use merge::merge_thread;
use output::output_thread;
use process::{analyze_thread, process_thread};

fn check_join(j: ScopedJoinHandle<anyhow::Result<()>>, s: &str) -> bool {
    if let Err(e) = j
        .join()
        .unwrap_or_else(|_| panic!("Error joining {s} thread"))
    {
        error!("{:?}", e);
        true
    } else {
        false
    }
}
fn merge_pipeline(cfg: Config) -> bool {
    let nt = cfg.threads();
    trace!("Running merge pipeline with {nt} threads");

    let mut error = false;

    thread::scope(|scope| {
        // Channel used to send files to read and merge thread
        let (sd, rx) = bounded(2);

        // Channel to send merged datasets for analysis
        let (sd_data, rx_data) = bounded(nt * 2);

        // Channel used to send results to output thread
        let (sd_res, rc_res) = unbounded();

        // Start output thread
        let cfg1 = &cfg;
        let output_task = scope.spawn(move |_| output_thread(cfg1, rc_res));

        // Add merge thread
        let cfg1 = &cfg;
        let merge_task = scope.spawn(move |_| merge_thread(cfg1, rx, sd_data));

        let mut process_tasks = Vec::with_capacity(nt);
        for ix in 0..nt {
            let rx1 = rx_data.clone();
            let sd_res1 = sd_res.clone();
            let cfg = &cfg;
            process_tasks.push(scope.spawn(move |_| analyze_thread(cfg, ix, rx1, sd_res1)));
        }

        drop(rx_data);
        drop(sd_res);

        for p in cfg.input_files() {
            sd.send(p)
                .expect("Error sending input file to process threads")
        }
        drop(sd);
        // Wait for merge thread
        error = check_join(merge_task, "merge thread");
        // ... and process threads
        for jh in process_tasks.drain(..) {
            error = error || check_join(jh, "process thread")
        }
        // ...and output thread
        error = error || check_join(output_task, "output thread")
    })
    .expect("Error in scope generation");

    error
}

fn std_pipeline(cfg: Config) -> bool {
    let nt = cfg.threads();
    trace!("Running standard pipeline with {nt} threads");
    let mut error = false;

    thread::scope(|scope| {
        // Channel used to send files to process threads
        let (sd, rx) = bounded(nt * 2);

        // Channel used to send results to output thread
        let (sd_res, rc_res) = unbounded();

        // Start output thread
        let cfg1 = &cfg;
        let output_task = scope.spawn(move |_| output_thread(cfg1, rc_res));

        let mut process_tasks = Vec::with_capacity(nt);
        for ix in 0..nt {
            let rx1 = rx.clone();
            let sd_res1 = sd_res.clone();
            let cfg = &cfg;
            process_tasks.push(scope.spawn(move |_| process_thread(cfg, ix, rx1, sd_res1)));
        }

        drop(rx);
        drop(sd_res);

        for p in cfg.input_files() {
            sd.send(p)
                .expect("Error sending input file to process threads")
        }
        drop(sd);
        // Wait for process threads
        for jh in process_tasks.drain(..) {
            error = error || check_join(jh, "process thread")
        }
        // ...and output thread
        error = error || check_join(output_task, "output thread")
    })
    .expect("Error in scope generation");

    error
}

fn main() -> anyhow::Result<()> {
    let cfg = cli::handle_cli()?;

    if if cfg.merge_key().is_none() {
        std_pipeline(cfg)
    } else {
        merge_pipeline(cfg)
    } {
        Err(anyhow!("Error occurred during processing"))
    } else {
        Ok(())
    }
}
