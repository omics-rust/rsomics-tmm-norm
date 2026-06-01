#!/usr/bin/env python3
"""Deterministic gene x sample count matrix for TMM tests/benches.

Negative-binomial-ish counts with per-sample depth scaling and a fraction of
zeros, so TMM's trim + zero-exclusion paths are exercised. Fixed seed -> stable.

Usage: mkmatrix.py n_genes n_samples seed > counts.tsv
"""
import sys
import random


def main() -> None:
    n_genes = int(sys.argv[1])
    n_samples = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 42
    rng = random.Random(seed)

    depth = [0.5 + 1.5 * (j / max(1, n_samples - 1)) for j in range(n_samples)]
    base = [rng.lognormvariate(3.0, 1.8) for _ in range(n_genes)]

    out = sys.stdout
    out.write("gene\t" + "\t".join(f"S{j}" for j in range(n_samples)) + "\n")
    for i in range(n_genes):
        cells = [f"G{i}"]
        mu = base[i]
        for j in range(n_samples):
            lam = mu * depth[j]
            # gamma-poisson mixture (NB) with dispersion 0.2
            shape = 5.0
            g = rng.gammavariate(shape, lam / shape) if lam > 0 else 0.0
            k = 0
            if g > 0:
                # poisson via simple inversion (g small enough here on average)
                l = pow(2.718281828459045, -g)
                p = 1.0
                while True:
                    p *= rng.random()
                    if p <= l:
                        break
                    k += 1
            cells.append(str(k))
        out.write("\t".join(cells) + "\n")


if __name__ == "__main__":
    main()
