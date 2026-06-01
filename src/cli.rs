use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Section};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-tmm-norm", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub counts: PathBuf,
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        let n = rsomics_tmm_norm::run(&self.counts, &mut out)?;
        if !self.common.quiet {
            eprintln!("{n} TMM normalization factors computed");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "TMM (trimmed mean of M-values) normalization factors for a count matrix.",
    origin: None,
    usage_lines: &["<counts.tsv> [-o factors.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[FlagSpec {
            short: Some('o'),
            long: "output",
            aliases: &[],
            value: Some("<path>"),
            type_hint: Some("String"),
            required: false,
            default: Some("-"),
            description: "Output TSV (sample<TAB>norm.factor); '-' is stdout.",
            why_default: Some("stream to stdout by default."),
        }],
    }],
    examples: &[Example {
        description: "Compute TMM factors for a gene x sample count matrix",
        command: "rsomics-tmm-norm counts.tsv -o factors.tsv",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
