#!/usr/bin/env python3
"""Discordance between the three tGBS libraries of the control accession
PI 490282, broken down by genotype class and by read depth.

usage: replicate_discordance.py ALLELE_COUNT_DIR BIM OUT_DIR

ALLELE_COUNT_DIR holds the provider's per-library allele count files
(<library>.allele_counts.txt: CHR POS REF ALT REF_COUNT ALT_COUNT OTHER_COUNT
TOTAL GENOTYPE, genotype 0 missing, 1 homozygous reference, 2 homozygous
alternative, 3 heterozygous). BIM is the plink .bim of the retained marker set;
the comparison is reported on all delivered sites and on the retained sites.
Standard library only, so it runs in any Python 3."""
import csv
import os
import sys

LIBS = ("PI_490282", "PI_490282C2", "PI_490282C3")
DEPTH_BINS = ((0, 4, "2-4"), (5, 9, "5-9"), (10, 19, "10-19"), (20, 10 ** 9, "20 or more"))


def read_lib(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            rows[(f[0], int(f[1]))] = (int(f[7]), int(f[8]))      # total depth, genotype
    return rows


def compare(a, b, keys, pair, site_set):
    both = hom_hom = het_hom = het_both = het_any = 0
    depth_conc = [0, 0]
    depth_disc = [0, 0]
    by_depth = {lab: [0, 0, 0, 0] for _, _, lab in DEPTH_BINS}       # both called, discordant, any het, both het
    for k in keys:
        (da, ga), (db, gb) = a[k], b[k]
        if ga == 0 or gb == 0:
            continue
        both += 1
        d = (da + db) / 2
        disc = ga != gb
        if disc:
            depth_disc[0] += d; depth_disc[1] += 1
            if ga in (1, 2) and gb in (1, 2):
                hom_hom += 1
            else:
                het_hom += 1
        else:
            depth_conc[0] += d; depth_conc[1] += 1
        if ga == 3 or gb == 3:
            het_any += 1
            if ga == 3 and gb == 3:
                het_both += 1
        for lo, hi, lab in DEPTH_BINS:
            if lo <= d <= hi:
                by_depth[lab][0] += 1
                by_depth[lab][1] += disc
                if ga == 3 or gb == 3:
                    by_depth[lab][2] += 1
                    by_depth[lab][3] += 1 if (ga == 3 and gb == 3) else 0
                break
    summary = dict(PAIR=pair, SITE_SET=site_set, N_BOTH_CALLED=both,
                   CONCORDANCE=round(1 - (hom_hom + het_hom) / both, 4),
                   HOM_HOM_DISCORDANT=hom_hom, HET_HOM_DISCORDANT=het_hom,
                   HET_CONCORDANCE=round(het_both / het_any, 4) if het_any else "",
                   MEAN_DEPTH_CONCORDANT=round(depth_conc[0] / depth_conc[1], 1),
                   MEAN_DEPTH_DISCORDANT=round(depth_disc[0] / depth_disc[1], 1))
    depth_rows = []
    for _, _, lab in DEPTH_BINS:
        n, nd, ha, hb = by_depth[lab]
        depth_rows.append(dict(PAIR=pair, SITE_SET=site_set, MEAN_DEPTH=lab, N_BOTH_CALLED=n,
                               CONCORDANCE=round(1 - nd / n, 4) if n else "",
                               HET_CONCORDANCE=round(hb / ha, 4) if ha else ""))
    return summary, depth_rows


def main():
    ac_dir, bim, out_dir = sys.argv[1:4]
    libs = {n: read_lib(os.path.join(ac_dir, n + ".allele_counts.txt")) for n in LIBS}
    all_keys = sorted(set.intersection(*(set(v) for v in libs.values())))
    retained = set()
    cmap = {}
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ref", "chr_map.txt")) as fh:
        for line in fh:
            a, b = line.split()
            cmap[b] = a
    with open(bim) as fh:
        for line in fh:
            f = line.split()
            retained.add((cmap.get(f[0], f[0]), int(f[3])))
    ret_keys = [k for k in all_keys if k in retained]
    summaries, depth_rows = [], []
    for i in range(3):
        for j in range(i + 1, 3):
            pair = LIBS[i] + " vs " + LIBS[j]
            for keys, label in ((all_keys, "all delivered sites"), (ret_keys, "retained sites")):
                s, d = compare(libs[LIBS[i]], libs[LIBS[j]], keys, pair, label)
                summaries.append(s); depth_rows.extend(d)
    os.makedirs(out_dir, exist_ok=True)
    for name, rows in (("review_replicate_discordance.csv", summaries),
                       ("review_replicate_concordance_by_depth.csv", depth_rows)):
        with open(os.path.join(out_dir, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
    for s in summaries:
        print(s)
    for d in depth_rows:
        print(d)


if __name__ == "__main__":
    main()
