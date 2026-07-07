use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-tmm-norm"))
}

fn golden(n: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), n)
}

fn parse(s: &str) -> Vec<(String, f64)> {
    s.trim()
        .lines()
        .skip(1)
        .map(|l| {
            let mut it = l.split('\t');
            let sample = it.next().unwrap().to_string();
            let f: f64 = it.next().unwrap().parse().unwrap();
            (sample, f)
        })
        .collect()
}

fn run_ours() -> Vec<(String, f64)> {
    let out = Command::new(ours())
        .arg(golden("counts.tsv"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    parse(&String::from_utf8(out.stdout).unwrap())
}

fn assert_close(mine: &[(String, f64)], theirs: &[(String, f64)]) {
    assert_eq!(mine.len(), theirs.len(), "sample count mismatch");
    let mut max_dev = 0.0f64;
    for (a, b) in mine.iter().zip(theirs.iter()) {
        assert_eq!(a.0, b.0, "sample id mismatch");
        let rel = (a.1 - b.1).abs() / b.1.abs().max(1e-12);
        max_dev = max_dev.max(rel);
        assert!(
            rel < 1e-6,
            "sample {}: ours={} oracle={} rel={rel:e}",
            a.0,
            a.1,
            b.1
        );
    }
    eprintln!("max relative deviation = {max_dev:e}");
}

// Always-on: diff against the committed edgeR golden so CI validates without R.
#[test]
fn matches_committed_golden() {
    let golden_out = std::fs::read_to_string(golden("factors.golden.tsv")).unwrap();
    assert_close(&run_ours(), &parse(&golden_out));
}

fn run_expecting_error(fixture: &str) -> String {
    let out = Command::new(ours()).arg(golden(fixture)).output().unwrap();
    assert!(
        !out.status.success(),
        "expected {fixture} to fail, but it succeeded"
    );
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(
        !stderr.contains("panicked"),
        "{fixture} panicked instead of failing loudly: {stderr}"
    );
    stderr
}

// edgeR errors ("missing value where TRUE/FALSE needed") on an all-zero library
// column. We must reject it loudly, not panic on a 0/0 NaN in the quantile sort.
#[test]
fn all_zero_column_errors_not_panic() {
    let stderr = run_expecting_error("allzero_column.tsv");
    assert!(
        stderr.contains("zero total count") || stderr.contains("missing value"),
        "unexpected error text: {stderr}"
    );
}

// edgeR errors ("NA counts not permitted") on a non-finite count literal. We must
// reject it at parse, not accept it and panic on a NaN in a downstream sort.
#[test]
fn non_finite_literal_errors_not_panic() {
    let stderr = run_expecting_error("nan_literal.tsv");
    assert!(
        stderr.contains("NA counts not permitted"),
        "unexpected error text: {stderr}"
    );
}

// Live differential vs edgeR via `conda run -n r-bioc Rscript`. Loud-skips when
// the r-bioc env is unavailable (e.g. CI runners with no Bioconductor).
#[test]
fn matches_edger_oracle() {
    let conda = match which_conda() {
        Some(c) => c,
        None => {
            eprintln!("SKIP matches_edger_oracle: no conda r-bioc env available");
            return;
        }
    };

    let oracle = format!("{}/tests/tmm_oracle.R", env!("CARGO_MANIFEST_DIR"));
    let ref_out = Command::new(&conda)
        .args([
            "run",
            "-n",
            "r-bioc",
            "Rscript",
            &oracle,
            &golden("counts.tsv"),
        ])
        .output()
        .unwrap();
    if !ref_out.status.success() {
        eprintln!(
            "SKIP matches_edger_oracle: oracle failed: {}",
            String::from_utf8_lossy(&ref_out.stderr)
        );
        return;
    }

    let theirs = parse(&String::from_utf8(ref_out.stdout).unwrap());
    assert_close(&run_ours(), &theirs);
}

fn which_conda() -> Option<String> {
    for c in [
        "conda".to_string(),
        format!(
            "{}/miniconda3/bin/conda",
            std::env::var("HOME").unwrap_or_default()
        ),
    ] {
        let ok = Command::new(&c)
            .args(["run", "-n", "r-bioc", "Rscript", "-e", "cat(1)"])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(c);
        }
    }
    None
}
