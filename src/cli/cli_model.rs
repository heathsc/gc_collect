use std::path::PathBuf;

use clap::{builder::PossibleValue, command, value_parser, Arg, ArgAction, Command, ValueEnum};

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
            Arg::new("merge")
                .short('m')
                .action(ArgAction::SetTrue)
                .long("merge")
                .help("Merge results in grouups"),
        )
        .arg(
            Arg::new("merge_by")
                .short('M')
                .long("merge-by")
                .value_name("MERGE KEY")
                .value_parser(value_parser!(MergeKey))
                .ignore_case(true)
                .help("Set merge key"),
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
            Arg::new("kmers")
                .long("kmers")
                .short('k')
                .value_parser(value_parser!(PathBuf))
                .value_name("KM FILE")
                .help("Input KM file with kmers for coverage estimation"),
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
            Arg::new("no_header")
                .short('H')
                .long("no-header")
                .action(ArgAction::SetTrue)
                .help("Do not add header to output file"),
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

#[derive(Debug, Clone, Copy)]
pub enum MergeKey {
    Default,
    Sample,
    Barcode,
    Library,
    Fli,
}

impl ValueEnum for MergeKey {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            Self::Default,
            Self::Sample,
            Self::Barcode,
            Self::Library,
            Self::Fli,
        ]
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            Self::Default => Some(PossibleValue::new("default")),
            Self::Sample => Some(PossibleValue::new("sample")),
            Self::Barcode => Some(PossibleValue::new("barcode")),
            Self::Library => Some(PossibleValue::new("library")),
            Self::Fli => Some(PossibleValue::new("fli")),
        }
    }
}
