#!/usr/bin/env python3
"""Read depth and allele balance at the sites removed by the heterozygous call
rate ceiling, against the retained sites.

Collapsed paralogs are read at roughly twice the depth of a single-copy site,
because the reads of two loci pile up at one position, and they are called
heterozygous in most accessions whatever their breeding history. This script
takes the per-library allele counts delivered with the genotypes (one file per
library: CHR, POS, REF, ALT, REF_COUNT, ALT_COUNT, OTHER_COUNT, TOTAL, GENOTYPE),
restricts them to the analysis accessions and the chromosome-anchored sites,
and compares the sites that the 0.55 ceiling removed with the sites it retained,
on depth per called genotype, depth relative to each library's own mean, depth
at heterozygous calls, and the share of reference reads at heterozygous calls.
The same quantities are given in bins of the heterozygous call rate.

usage: site_depth_by_filter.py ALLELE_COUNT_DIR GENOTYPE_TXT KEEP_SAMPLES OUT_DIR
"""
import os
import sys

import numpy as np
import pandas as pd

CHROM_PREFIX = "NC_0354"            # Ha1 to Ha17 of HanXRQr2.0 (NC_035433.2 to NC_035449.2)
MAX_SITE_MISS, MIN_MAF, MAX_SITE_HET = 0.50, 0.05, 0.55
BINS = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.55), (0.55, 1.01)]


def main(acdir, geno_txt, keep_file, outdir):
    keep = [l.split()[-1] for l in open(keep_file) if l.strip()]
    sites = pd.read_csv(geno_txt, sep="\t", comment="#", usecols=["Chr", "Pos"], dtype={"Chr": str, "Pos": np.int64})
    key = pd.MultiIndex.from_frame(sites)
    n_sites = len(sites)
    called = np.zeros(n_sites, dtype=np.int32); het = np.zeros(n_sites, dtype=np.int32)
    alt_alleles = np.zeros(n_sites, dtype=np.int64)
    depth_called = np.zeros(n_sites); depth_rel = np.zeros(n_sites)
    depth_het = np.zeros(n_sites); ref_share_het = np.zeros(n_sites)
    for s in keep:
        f = os.path.join(acdir, s + ".allele_counts.txt")
        d = pd.read_csv(f, sep="\t", comment="#", header=None,
                        names=["CHR", "POS", "REF", "ALT", "REF_COUNT", "ALT_COUNT", "OTHER", "TOTAL", "GENOTYPE"],
                        usecols=["CHR", "POS", "REF_COUNT", "ALT_COUNT", "TOTAL", "GENOTYPE"], dtype={"CHR": str})
        d = d.set_index(["CHR", "POS"]).reindex(key)
        g = d["GENOTYPE"].fillna(0).to_numpy().astype(np.int8)
        tot = d["TOTAL"].fillna(0).to_numpy().astype(float)
        rc = d["REF_COUNT"].fillna(0).to_numpy().astype(float); ac = d["ALT_COUNT"].fillna(0).to_numpy().astype(float)
        c = g != 0; h = g == 3
        mean_depth = tot[c].mean()
        called += c; het += h
        alt_alleles += np.where(g == 2, 2, 0) + np.where(h, 1, 0)
        depth_called += np.where(c, tot, 0.0); depth_rel += np.where(c, tot / mean_depth, 0.0)
        depth_het += np.where(h, tot, 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            ref_share_het += np.where(h & (rc + ac > 0), rc / np.maximum(rc + ac, 1), 0.0)
        print(s, "mean depth %.1f" % mean_depth, file=sys.stderr)
    n = len(keep)
    chrom_ok = sites["Chr"].str.startswith(CHROM_PREFIX).to_numpy()
    call_rate = called / n
    p = alt_alleles / np.maximum(2 * called, 1); maf = np.minimum(p, 1 - p)
    het_rate = het / np.maximum(called, 1)
    base = chrom_ok & (call_rate >= 1 - MAX_SITE_MISS) & (maf >= MIN_MAF)
    removed = base & (het_rate > MAX_SITE_HET); kept = base & (het_rate <= MAX_SITE_HET)
    site_depth = depth_called / np.maximum(called, 1)

    def row(label, m):
        return dict(SITE_CLASS=label, N_SITES=int(m.sum()),
                    MEAN_DEPTH_PER_CALLED_GENOTYPE=round(depth_called[m].sum() / called[m].sum(), 2),
                    MEDIAN_SITE_MEAN_DEPTH=round(float(np.median(site_depth[m])), 2),
                    MEAN_RELATIVE_DEPTH=round(depth_rel[m].sum() / called[m].sum(), 3),
                    MEAN_DEPTH_AT_HETEROZYGOUS_CALLS=round(depth_het[m].sum() / max(het[m].sum(), 1), 2),
                    REF_SHARE_AT_HETEROZYGOUS_CALLS=round(ref_share_het[m].sum() / max(het[m].sum(), 1), 3),
                    MEAN_HET_CALL_RATE=round(float(het_rate[m].mean()), 3))
    rows = [row("retained (heterozygous call rate <= 0.55)", kept),
            row("removed by the 0.55 ceiling", removed)]
    for lo, hi in BINS:
        m = base & (het_rate > lo) & (het_rate <= hi) if lo > 0 else base & (het_rate <= hi)
        rows.append(row("heterozygous call rate %s to %s" % (lo, min(hi, 1.0)), m))
    out = pd.DataFrame(rows)
    out.loc[len(out)] = dict(SITE_CLASS="ratio removed / retained", N_SITES="",
                             **{c: round(out.loc[1, c] / out.loc[0, c], 3) for c in out.columns[2:]})
    os.makedirs(outdir, exist_ok=True)
    out.to_csv(os.path.join(outdir, "review_site_depth_by_filter.csv"), index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
