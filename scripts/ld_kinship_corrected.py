#!/usr/bin/env python3
"""Linkage disequilibrium corrected for kinship (r2_V) and for kinship and
structure (r2_VS), after Mangin et al. (2012), against the uncorrected r2 on
the same marker pairs.

For every pair of retained markers within 1 Mb on the same chromosome the
squared correlation is computed three ways on the imputed genotypes:
  r2     on the centered genotype vectors (the estimator of the paper);
  r2_V   on the vectors transformed by V^-1/2, where V is the genomic
         relationship matrix of the accessions in the set, with the mean
         estimated by generalized least squares (kinship correction);
  r2_VS  as r2_V after regressing the transformed vectors on the transformed
         ancestry proportions (kinship and structure correction; whole panel).
Means are taken in 2-kb bins of distance, as in the empirical curves of the
paper, for the whole panel and for each cluster with at least 20 accessions,
with the kinship matrix restricted to the accessions of the set.

usage: ld_kinship_corrected.py IMPUTED_GENOTYPES KINSHIP_CSV ACCESSION_TABLE Q_TABLE RETAINED_SITES OUT_DIR
  IMPUTED_GENOTYPES  the provider's imputed matrix (Chr, Pos, REF, ALT, one column per library)
  KINSHIP_CSV        the genomic relationship matrix (Data file 1)
  ACCESSION_TABLE    Table S1 (SAMPLE, CLUSTER)
  Q_TABLE            Table S6 (SAMPLE, Q1..Q5)
  RETAINED_SITES     two columns, chromosome and position, of the retained markers
"""
import os
import sys

import numpy as np
import pandas as pd

BIN = 2000
MAXDIST = 1_000_000
NBIN = MAXDIST // BIN
AT = {"R2_1KB": 1, "R2_5KB": 5, "R2_10KB": 10, "R2_20KB": 20, "R2_50KB": 50, "R2_100KB": 100, "R2_200KB": 200, "R2_500KB": 500}
MIN_CLUSTER = 20
MIN_MAF = 0.05


def whitener(V, floor=1e-4):
    """V^-1/2 by eigendecomposition, with small eigenvalues raised to a floor."""
    w, U = np.linalg.eigh(V)
    w = np.maximum(w, floor * w.max())
    return (U * (1.0 / np.sqrt(w))) @ U.T


def residualize(Z, C):
    """Regress the columns of Z on the covariate columns C."""
    beta, *_ = np.linalg.lstsq(C, Z, rcond=None)
    return Z - C @ beta


def standardize(Z):
    Z = Z - Z.mean(0)
    sd = Z.std(0)
    sd[sd == 0] = np.inf
    return (Z / sd).astype(np.float32)


def curve(Z, chrom, pos):
    """Sum and count of r2 in distance bins, over pairs within MAXDIST on each chromosome.
    Z is samples x markers, standardized (mean 0, variance 1 per column)."""
    n = Z.shape[0]
    s = np.zeros(NBIN); c = np.zeros(NBIN, dtype=np.int64)
    for ch in np.unique(chrom):
        idx = np.where(chrom == ch)[0]
        idx = idx[np.argsort(pos[idx], kind="stable")]
        p = pos[idx].astype(np.int64); M = Z[:, idx]
        B = 1024
        for a in range(0, len(idx), B):
            b = min(a + B, len(idx))
            e = np.searchsorted(p, p[b - 1] + MAXDIST, side="right")
            R = (M[:, a:b].T @ M[:, a:e]) / n                 # correlations, block x window
            d = p[a:e][None, :] - p[a:b][:, None]              # distance, positive for pairs ahead
            ok = (d > 0) & (d < MAXDIST)
            k = (d[ok] // BIN).astype(np.int64)
            s += np.bincount(k, weights=(R[ok] ** 2).astype(np.float64), minlength=NBIN)
            c += np.bincount(k, minlength=NBIN)
    return s, c


def summarize(label, estimator, s, c):
    r2 = s / np.maximum(c, 1)
    row = dict(SET=label, ESTIMATOR=estimator, PAIRS=int(c.sum()))
    for col, d in AT.items():
        row[col] = round(float(r2[d * 1000 // BIN]), 4)
    bg = slice(800_000 // BIN, NBIN)
    row["R2_800KB_TO_1MB"] = round(float(s[bg].sum() / max(c[bg].sum(), 1)), 4)
    return row, r2


def main(geno, kin, acc_tab, q_tab, sites, outdir):
    acc = pd.read_csv(acc_tab)
    K = pd.read_csv(kin, index_col=0)
    Q = pd.read_csv(q_tab).set_index("SAMPLE")[["Q1", "Q2", "Q3", "Q4", "Q5"]]
    samples = list(acc["SAMPLE"])
    keep = pd.read_csv(sites, sep=r"\s+", header=None, names=["CHR", "POS"], dtype={"CHR": str})
    keyset = set(zip(keep["CHR"], keep["POS"]))
    g = pd.read_csv(geno, sep="\t", comment="#", dtype={"Chr": str})
    sel = [i for i, (c, p) in enumerate(zip(g["Chr"], g["Pos"])) if (c, p) in keyset]
    g = g.iloc[sel].reset_index(drop=True)
    chrom = g["Chr"].to_numpy(); pos = g["Pos"].to_numpy()
    X = g[samples].to_numpy(dtype=np.int8).T                    # samples x markers, codes 1 ref, 2 alt, 3 het
    X = np.where(X == 1, 0, np.where(X == 2, 2, np.where(X == 3, 1, -1))).astype(np.float32)
    assert (X >= 0).all(), "missing genotypes in the imputed matrix"
    print("%d accessions x %d retained markers" % X.shape, file=sys.stderr)

    sets = [("All accessions", samples)]
    for cl, grp in acc.groupby("CLUSTER"):
        if cl != "admixed" and len(grp) >= MIN_CLUSTER:
            sets.append((cl, list(grp["SAMPLE"])))
    rows, binned = [], []
    for label, ids in sets:
        ii = [samples.index(s) for s in ids]
        Xs = X[ii]
        p = Xs.mean(0) / 2                                      # markers below MIN_MAF within the set are dropped,
        poly = np.minimum(p, 1 - p) >= MIN_MAF                  # as PopLDdecay did with -MAF for each subpopulation
        Xs = Xs[:, poly]; ch_s = chrom[poly]; pos_s = pos[poly]
        print("%s: %d accessions, %d markers with MAF >= %.2f" % (label, len(ids), poly.sum(), MIN_MAF), file=sys.stderr)
        V = K.loc[ids, ids].to_numpy(dtype=np.float64)
        T = whitener(V)
        one = np.ones((len(ids), 1))
        est = {"r2": standardize(Xs),
               "r2_V": standardize(residualize(T @ Xs, T @ one))}
        if label == "All accessions":
            S = Q.loc[ids].to_numpy()[:, :-1]                   # four ancestry columns, the fifth is implied
            est["r2_VS"] = standardize(residualize(T @ Xs, np.hstack([T @ one, T @ S])))
        for name, Z in est.items():
            s, c = curve(Z, ch_s, pos_s)
            row, r2 = summarize(label, name, s, c)
            rows.append(row)
            binned.extend(dict(SET=label, ESTIMATOR=name, DISTANCE_KB=(k + 0.5) * BIN / 1000, R2=round(float(r2[k]), 5), PAIRS=int(c[k]))
                          for k in range(NBIN) if c[k] > 0)
            print(label, name, {k: row[k] for k in ("R2_5KB", "R2_100KB", "R2_800KB_TO_1MB")}, file=sys.stderr)
    out = pd.DataFrame(rows)
    ref = out[out["ESTIMATOR"] == "r2"].set_index("SET")
    for col in list(AT) + ["R2_800KB_TO_1MB"]:
        out["RATIO_" + col] = (out[col] / out["SET"].map(ref[col])).round(3)
    os.makedirs(outdir, exist_ok=True)
    out.to_csv(os.path.join(outdir, "review_ld_kinship_corrected.csv"), index=False)
    pd.DataFrame(binned).to_csv(os.path.join(outdir, "review_ld_kinship_corrected_binned.csv"), index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:7]))
