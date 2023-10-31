use std::path::PathBuf;

use clap::{command, value_parser, Arg, ArgAction, Command};

use crate::utils::LogLevel;

pub(super) fn cli_model() -> Command {
    command!()
        .arg(
            Arg::new("timestamp")
                .short('X')
                .long("timestamp")
                .value_parser(value_parser!(stderrlog::Timestamp))
                .value_name("GRANULARITY")
                .default_value("none")
                .help("Prepend log entries with a timestamp"),
        )
        .arg(
            Arg::new("loglevel")
                .short('l')
                .long("loglevel")
                .value_name("LOGLEVEL")
                .value_parser(value_parser!(LogLevel))
                .ignore_case(true)
                .default_value("info")
                .help("Set log level"),
        )
        .arg(
            Arg::new("quiet")
                .action(ArgAction::SetTrue)
                .long("quiet")
                .conflicts_with("loglevel")
                .help("Silence all output"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_parser(value_parser!(u64).range(1..))
                .value_name("INT")
                .help("Set number of process threads [default: number of available cores]"),
        )
        .arg(
            Arg::new("ref")
                .short('r')
                .long("reference-json")
                .value_parser(value_parser!(PathBuf))
                .value_name("FILE")
                .help("Reference JSON file produced by analyze_ref_gc"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_parser(value_parser!(PathBuf))
                .value_name("OUTPUT")
                .help("Main output file [default: <stdout>]"),
        )
        .arg(
            Arg::new("input")
                .value_parser(value_parser!(PathBuf))
                .value_name("INPUT")
                .num_args(1..)
                .required(true)
                .help("Input JSON file(s) from fastq_gc"),
        )
}
