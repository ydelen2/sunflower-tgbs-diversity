#!/usr/bin/env python3
"""
Per-accession QC table from a *.snps.genotype.txt matrix.

Columns: call rate, missing rate, hom-REF / hom-ALT / het counts, Ho over
called sites, and the same restricted to sites called in >= --core-call of
the retained accessions (CORE_HO).

Usage
  python3 sample_stats.py -i Sunflower.MCR50.snps.genotype.txt \
      -o results/01_qc/sample_stats.tsv --core-call 0.90
"""
import argparse, sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--core-call", type=float, default=0.90,
                    help="site call-rate threshold defining the common site subset")
    ap.add_argument("--chr-map", help="restrict to sequences listed here (chromosomes only)")
    ap.add_argument("--keep", help="file with one sample ID per line; restrict to these accessions")
    a = ap.parse_args()

    keep_chr = None
    if a.chr_map:
        keep_chr = {l.split()[0] for l in open(a.chr_map) if l.strip() and not l.startswith("#")}

    fh = open(a.input)
    line = fh.readline()
    while line.startswith("#"):
        line = fh.readline()
    all_samples = line.rstrip("\n").split("\t")[4:]
    if a.keep:
        want = {l.strip() for l in open(a.keep) if l.strip()}
        col = [i for i, s_ in enumerate(all_samples) if s_ in want]
    else:
        col = list(range(len(all_samples)))
    samples = [all_samples[i] for i in col]
    n = len(samples)

    miss = [0] * n; hom_r = [0] * n; hom_a = [0] * n; het = [0] * n
    c_miss = [0] * n; c_hom = [0] * n; c_het = [0] * n
    n_site = n_core = 0

    for line in fh:
        f = line.rstrip("\n").split("\t")
        if keep_chr is not None and f[0] not in keep_chr:
            continue
        gall = f[4:]
        g = [gall[i] for i in col]
        n_site += 1
        nm = g.count("0")
        core = (n - nm) / n >= a.core_call
        if core:
            n_core += 1
        for i, c in enumerate(g):
            if c == "0":
                miss[i] += 1
                if core: c_miss[i] += 1
            elif c == "3":
                het[i] += 1
                if core: c_het[i] += 1
            elif c == "1":
                hom_r[i] += 1
                if core: c_hom[i] += 1
            else:
                hom_a[i] += 1
                if core: c_hom[i] += 1

    with open(a.output, "w") as out:
        out.write("SAMPLE\tN_SITES\tN_MISS\tMISS_RATE\tCALL_RATE\tN_HOM_REF\tN_HOM_ALT"
                  "\tN_HET\tHO\tCORE_SITES\tCORE_MISS_RATE\tCORE_HO\n")
        for i, s in enumerate(samples):
            called = n_site - miss[i]
            ho = het[i] / called if called else float("nan")
            c_called = n_core - c_miss[i]
            c_ho = c_het[i] / c_called if c_called else float("nan")
            out.write(f"{s}\t{n_site}\t{miss[i]}\t{miss[i]/n_site:.4f}"
                      f"\t{called/n_site:.4f}\t{hom_r[i]}\t{hom_a[i]}\t{het[i]}\t{ho:.4f}"
                      f"\t{n_core}\t{c_miss[i]/n_core if n_core else float('nan'):.4f}"
                      f"\t{c_ho:.4f}\n")

    print(f"[sample_stats] {n} accessions, {n_site} sites, "
          f"{n_core} sites called in >= {a.core_call:.0%} of accessions", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
