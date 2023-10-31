#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

use crossbeam_channel::{bounded, unbounded};
use crossbeam_utils::thread;
use crossbeam_utils::thread::ScopedJoinHandle;

mod betabin;
mod cli;
mod gauss_legendre;
mod output;
mod process;
mod read;
mod reference;
mod simple_regression;
mod utils;

use output::output_thread;
use process::process_thread;

fn main() -> anyhow::Result<()> {
    let cfg = cli::handle_cli()?;
    let nt = cfg.threads();

    let mut error = false;

    let mut jn = |j: ScopedJoinHandle<anyhow::Result<()>>, s: &str| {
        if let Err(e) = j
            .join()
            .unwrap_or_else(|_| panic!("Error joining {s} thread"))
        {
            error!("{:?}", e);
            error = true
        }
    };

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
            jn(jh, "process thread");
        }
        // ...and output thread
        jn(output_task, "output thread")
    })
    .expect("Error in scope generation");

    if error {
        Err(anyhow!("Error occurred during processing"))
    } else {
        Ok(())
    }
}
