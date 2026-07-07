use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

const LOGRATIO_TRIM: f64 = 0.3;
const SUM_TRIM: f64 = 0.05;
const ACUTOFF: f64 = -1e10;

pub struct Matrix {
    pub samples: Vec<String>,
    pub counts: Vec<Vec<f64>>,
}

pub fn read_matrix(path: &Path) -> Result<Matrix> {
    let bytes = fs::read(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut lines = bytes.split(|&b| b == b'\n');

    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty count matrix".into()))?;
    let samples: Vec<String> = header
        .split(|&b| b == b'\t')
        .skip(1)
        .map(|s| String::from_utf8_lossy(s).into_owned())
        .collect();
    if samples.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "count matrix has no sample columns".into(),
        ));
    }

    let mut counts: Vec<Vec<f64>> = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let mut row = Vec::with_capacity(samples.len());
        for cell in line.split(|&b| b == b'\t').skip(1) {
            row.push(parse_count(cell)?);
        }
        if row.len() != samples.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "row has {} values but header has {} samples",
                row.len(),
                samples.len()
            )));
        }
        counts.push(row);
    }
    Ok(Matrix { samples, counts })
}

/// Counts are typically non-negative integers; take a no-allocation byte path
/// for that case and fall back to a full f64 parse for decimals.
fn parse_count(cell: &[u8]) -> Result<f64> {
    if cell.is_empty() {
        return Err(RsomicsError::InvalidInput("empty count cell".into()));
    }
    let mut acc: u64 = 0;
    let mut all_digits = true;
    for &b in cell {
        if b.is_ascii_digit() {
            acc = acc.wrapping_mul(10).wrapping_add((b - b'0') as u64);
        } else {
            all_digits = false;
            break;
        }
    }
    if all_digits {
        return Ok(acc as f64);
    }
    let s = std::str::from_utf8(cell)
        .map_err(|_| RsomicsError::InvalidInput("non-UTF8 count cell".into()))?;
    let v: f64 = s
        .parse()
        .map_err(|_| RsomicsError::InvalidInput(format!("non-numeric count '{s}'")))?;
    if !v.is_finite() {
        return Err(RsomicsError::InvalidInput("NA counts not permitted".into()));
    }
    if v < 0.0 {
        return Err(RsomicsError::InvalidInput(
            "negative counts not allowed".into(),
        ));
    }
    Ok(v)
}

/// edgeR calcNormFactors(method="TMM"). Returns one factor per sample,
/// scaled so their geometric mean is 1.
pub fn tmm_factors(m: &Matrix) -> Result<Vec<f64>> {
    let n_samples = m.samples.len();
    if n_samples == 0 {
        return Ok(Vec::new());
    }

    let kept: Vec<&Vec<f64>> = m
        .counts
        .iter()
        .filter(|row| row.iter().any(|&c| c > 0.0))
        .collect();

    if kept.is_empty() || n_samples == 1 {
        return Ok(vec![1.0; n_samples]);
    }

    let n_genes = kept.len();
    let mut lib_size = vec![0.0f64; n_samples];
    for row in &kept {
        for (j, &c) in row.iter().enumerate() {
            lib_size[j] += c;
        }
    }

    // An all-zero library column makes the 0.75-quantile / lib.size ratio 0/0 = NaN,
    // which edgeR propagates into a "missing value where TRUE/FALSE needed" error.
    if let Some(j) = lib_size.iter().position(|&s| s <= 0.0) {
        return Err(RsomicsError::InvalidInput(format!(
            "sample '{}' has zero total count: missing value where TRUE/FALSE needed",
            m.samples[j]
        )));
    }

    let ref_col = reference_column(&kept, &lib_size, n_genes, n_samples);

    let mut obs = vec![0.0f64; n_genes];
    let reference: Vec<f64> = (0..n_genes).map(|g| kept[g][ref_col]).collect();

    let mut f = vec![0.0f64; n_samples];
    for (j, fj) in f.iter_mut().enumerate() {
        for (g, o) in obs.iter_mut().enumerate() {
            *o = kept[g][j];
        }
        *fj = calc_factor_tmm(&obs, &reference, lib_size[j], lib_size[ref_col]);
    }

    let log_mean: f64 = f.iter().map(|x| x.ln()).sum::<f64>() / n_samples as f64;
    let scale = log_mean.exp();
    for fj in &mut f {
        *fj /= scale;
    }
    Ok(f)
}

fn reference_column(
    kept: &[&Vec<f64>],
    lib_size: &[f64],
    n_genes: usize,
    n_samples: usize,
) -> usize {
    let mut f75 = vec![0.0f64; n_samples];
    let mut col = vec![0.0f64; n_genes];
    for (j, slot) in f75.iter_mut().enumerate() {
        for (g, c) in col.iter_mut().enumerate() {
            *c = kept[g][j];
        }
        *slot = quantile_type7(&mut col, 0.75) / lib_size[j];
    }

    let mut sorted = f75.clone();
    let median = median_sorted(&mut sorted);
    if median < 1e-20 {
        // degenerate libraries: largest sqrt-mass column
        let mut sqrt_mass = vec![0.0f64; n_samples];
        for row in kept {
            for (j, m) in sqrt_mass.iter_mut().enumerate() {
                *m += row[j].sqrt();
            }
        }
        return argmax(&sqrt_mass);
    }

    let mean: f64 = f75.iter().sum::<f64>() / n_samples as f64;
    let mut best = 0usize;
    let mut best_dist = f64::INFINITY;
    for (j, &v) in f75.iter().enumerate() {
        let d = (v - mean).abs();
        if d < best_dist {
            best_dist = d;
            best = j;
        }
    }
    best
}

fn argmax(x: &[f64]) -> usize {
    x.iter()
        .enumerate()
        .fold((0usize, f64::NEG_INFINITY), |(bi, bm), (i, &v)| {
            if v > bm { (i, v) } else { (bi, bm) }
        })
        .0
}

fn calc_factor_tmm(obs: &[f64], reference: &[f64], n_o: f64, n_r: f64) -> f64 {
    let mut log_r = Vec::with_capacity(obs.len());
    let mut abs_e = Vec::with_capacity(obs.len());
    let mut var = Vec::with_capacity(obs.len());

    for (&o, &r) in obs.iter().zip(reference.iter()) {
        let lr = ((o / n_o) / (r / n_r)).log2();
        let ae = ((o / n_o).log2() + (r / n_r).log2()) / 2.0;
        let v = (n_o - o) / n_o / o + (n_r - r) / n_r / r;
        if lr.is_finite() && ae.is_finite() && ae > ACUTOFF {
            log_r.push(lr);
            abs_e.push(ae);
            var.push(v);
        }
    }

    if log_r.is_empty() {
        return 1.0;
    }
    if log_r.iter().fold(0.0f64, |m, &x| m.max(x.abs())) < 1e-6 {
        return 1.0;
    }

    let n = log_r.len();
    let lo_l = (n as f64 * LOGRATIO_TRIM).floor() + 1.0;
    let hi_l = n as f64 + 1.0 - lo_l;
    let lo_s = (n as f64 * SUM_TRIM).floor() + 1.0;
    let hi_s = n as f64 + 1.0 - lo_s;

    let rank_r = average_rank(&log_r);
    let rank_e = average_rank(&abs_e);

    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for i in 0..n {
        if rank_r[i] >= lo_l && rank_r[i] <= hi_l && rank_e[i] >= lo_s && rank_e[i] <= hi_s {
            let w = 1.0 / var[i];
            if w.is_finite() && (log_r[i] / var[i]).is_finite() {
                num += log_r[i] / var[i];
                den += w;
            }
        }
    }

    let f = if den == 0.0 { 0.0 } else { num / den };
    let f = if f.is_nan() { 0.0 } else { f };
    2.0f64.powf(f)
}

/// R's default `rank` with ties.method = "average".
fn average_rank(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| x[a].total_cmp(&x[b]));

    let mut ranks = vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && x[order[j]] == x[order[i]] {
            j += 1;
        }
        // ranks i..j (0-based) are tied; average of 1-based positions
        let avg = ((i + 1 + j) as f64) / 2.0;
        for &idx in &order[i..j] {
            ranks[idx] = avg;
        }
        i = j;
    }
    ranks
}

/// R quantile type 7 (the default): linear interpolation, h = (n-1)*p.
fn quantile_type7(x: &mut [f64], p: f64) -> f64 {
    let n = x.len();
    if n == 0 {
        return 0.0;
    }
    x.sort_by(|a, b| a.total_cmp(b));
    let h = (n - 1) as f64 * p;
    let lo = h.floor() as usize;
    let frac = h - lo as f64;
    if lo + 1 < n {
        x[lo] + frac * (x[lo + 1] - x[lo])
    } else {
        x[lo]
    }
}

fn median_sorted(x: &mut [f64]) -> f64 {
    x.sort_by(|a, b| a.total_cmp(b));
    let n = x.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 1 {
        x[n / 2]
    } else {
        (x[n / 2 - 1] + x[n / 2]) / 2.0
    }
}

pub fn write_factors(samples: &[String], factors: &[f64], output: &mut dyn Write) -> Result<()> {
    let mut out = BufWriter::new(output);
    writeln!(out, "sample\tnorm.factor").map_err(RsomicsError::Io)?;
    for (s, f) in samples.iter().zip(factors.iter()) {
        writeln!(out, "{s}\t{f:.10}").map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

pub fn run(counts_path: &Path, output: &mut dyn Write) -> Result<usize> {
    let m = read_matrix(counts_path)?;
    let factors = tmm_factors(&m)?;
    write_factors(&m.samples, &factors, output)?;
    Ok(m.samples.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn average_rank_handles_ties() {
        // R: rank(c(3,1,1,2)) == c(4.0, 1.5, 1.5, 3.0)
        let r = average_rank(&[3.0, 1.0, 1.0, 2.0]);
        assert_eq!(r, vec![4.0, 1.5, 1.5, 3.0]);
    }

    #[test]
    fn quantile_type7_matches_r() {
        // R: quantile(c(1,2,3,4), 0.75) == 3.25
        let mut v = vec![1.0, 2.0, 3.0, 4.0];
        assert!((quantile_type7(&mut v, 0.75) - 3.25).abs() < 1e-12);
        // R: quantile(c(10,20,30), 0.75) == 25
        let mut v = vec![30.0, 10.0, 20.0];
        assert!((quantile_type7(&mut v, 0.75) - 25.0).abs() < 1e-12);
    }

    #[test]
    fn identical_columns_give_unit_factors() {
        let m = Matrix {
            samples: vec!["a".into(), "b".into(), "c".into()],
            counts: vec![
                vec![100.0, 100.0, 100.0],
                vec![50.0, 50.0, 50.0],
                vec![10.0, 10.0, 10.0],
                vec![200.0, 200.0, 200.0],
            ],
        };
        for f in tmm_factors(&m).unwrap() {
            assert!((f - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn factors_have_unit_geometric_mean() {
        let m = Matrix {
            samples: vec!["a".into(), "b".into()],
            counts: vec![
                vec![100.0, 200.0],
                vec![50.0, 80.0],
                vec![10.0, 5.0],
                vec![300.0, 600.0],
                vec![5.0, 0.0],
            ],
        };
        let f = tmm_factors(&m).unwrap();
        let log_mean: f64 = f.iter().map(|x| x.ln()).sum::<f64>() / f.len() as f64;
        assert!(log_mean.abs() < 1e-9);
    }

    #[test]
    fn all_zero_column_errors_cleanly() {
        // edgeR: an all-zero library column -> "missing value where TRUE/FALSE needed".
        let m = Matrix {
            samples: vec!["S0".into(), "S1".into(), "S2".into()],
            counts: vec![
                vec![10.0, 0.0, 5.0],
                vec![20.0, 0.0, 8.0],
                vec![5.0, 0.0, 3.0],
            ],
        };
        let err = tmm_factors(&m).unwrap_err();
        assert!(matches!(err, RsomicsError::InvalidInput(_)));
    }

    #[test]
    fn non_finite_count_literal_rejected() {
        // edgeR: an NA/non-finite count -> "NA counts not permitted".
        assert!(parse_count(b"NaN").is_err());
        assert!(parse_count(b"inf").is_err());
        assert!(parse_count(b"-inf").is_err());
    }
}
