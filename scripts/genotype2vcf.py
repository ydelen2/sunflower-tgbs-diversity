#!/usr/bin/env python3
"""
Convert a *.snps.genotype.txt matrix into a minimal GT-only VCF.

The delivered *.snps.vcf declares ##INFO=<ID=AF,Number=A,...> but writes three
values per record, and uses AC as a FORMAT key; bcftools and GATK reject both.

Input  columns : Chr  Pos  REF  ALT  <sample1> <sample2> ...
Genotype codes : 0 = missing, 1 = hom REF, 2 = hom ALT, 3 = heterozygous

Usage
  python3 genotype2vcf.py -i Sunflower.MCR50.snps.genotype.txt \
      --chr-map ref/chr_map.txt --chrom-only \
      --max-site-miss 0.20 --min-maf 0.05 --max-site-het 0.20 \
      --keep keep_samples.txt --stats site_stats.tsv | bgzip -c > out.vcf.gz
"""
import argparse, sys, time

CODE = {"0": "./.", "1": "0/0", "2": "1/1", "3": "0/1"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", default="-")
    ap.add_argument("--chr-map", help="two-column TSV: old_name<TAB>new_name")
    ap.add_argument("--chrom-only", action="store_true",
                    help="keep only sequences listed in --chr-map (drops scaffolds/organelles)")
    ap.add_argument("--keep", help="file with one sample ID per line to retain")
    ap.add_argument("--max-site-miss", type=float, default=1.0)
    ap.add_argument("--min-maf", type=float, default=0.0)
    ap.add_argument("--max-site-het", type=float, default=1.0)
    ap.add_argument("--min-called", type=int, default=2,
                    help="minimum number of non-missing genotypes at a site")
    ap.add_argument("--stats", help="write per-site MAF/He/Ho/missing table here")
    ap.add_argument("--reference", default="HanXRQr2.0-SUNRISE")
    a = ap.parse_args()

    cmap = {}
    if a.chr_map:
        with open(a.chr_map) as fh:
            for ln in fh:
                if ln.strip() and not ln.startswith("#"):
                    old, new = ln.split()[:2]
                    cmap[old] = new

    keep = None
    if a.keep:
        with open(a.keep) as fh:
            keep = {ln.strip() for ln in fh if ln.strip()}

    fin = open(a.input)
    fout = sys.stdout if a.output == "-" else open(a.output, "w")
    fstat = open(a.stats, "w") if a.stats else None
    if fstat:
        fstat.write("CHROM\tPOS\tN_CALLED\tMISS\tMAF\tHE\tHO\tFIS\tPASS\tREASON\n")

    # ---- header ------------------------------------------------------------
    line = fin.readline()
    while line.startswith("#"):                       # skip the code-legend line
        line = fin.readline()
    hdr = line.rstrip("\n").split("\t")
    samples = hdr[4:]
    idx = [i for i, s in enumerate(samples)] if keep is None else \
          [i for i, s in enumerate(samples) if s in keep]
    kept_samples = [samples[i] for i in idx]
    if keep and len(kept_samples) != len(keep):
        missing = keep - set(kept_samples)
        print(f"[warn] {len(missing)} IDs in --keep not found in matrix: "
              f"{sorted(missing)[:5]}...", file=sys.stderr)

    fout.write("##fileformat=VCFv4.2\n")
    fout.write(f"##fileDate={time.strftime('%Y%m%d')}\n")
    fout.write("##source=genotype2vcf.py\n")
    fout.write(f"##reference={a.reference}\n")
    for old, new in cmap.items():
        fout.write(f'##contig=<ID={new}>\n')
    fout.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">\n')
    fout.write('##INFO=<ID=AF,Number=A,Type=Float,Description="ALT allele frequency">\n')
    fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
               + "\t".join(kept_samples) + "\n")

    n_in = n_out = 0
    drop_scaf = drop_miss = drop_maf = drop_het = 0

    for line in fin:
        f = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = f[0], f[1], f[2], f[3]
        n_in += 1

        if cmap:
            if chrom in cmap:
                chrom = cmap[chrom]
            elif a.chrom_only:
                drop_scaf += 1
                continue

        g = f[4:]
        codes = [g[i] for i in idx]
        n0 = n1 = n2 = n3 = 0
        for c in codes:
            if c == "1":   n1 += 1
            elif c == "2": n2 += 1
            elif c == "3": n3 += 1
            else:          n0 += 1
        n_called = n1 + n2 + n3
        n_tot = n_called + n0
        if n_called < a.min_called:
            drop_miss += 1
            continue

        miss = n0 / n_tot
        het = n3 / n_called
        p = (2 * n1 + n3) / (2 * n_called)          # REF allele freq
        q = 1.0 - p
        maf = p if p < q else q
        he = 2 * p * q
        fis = (1 - het / he) if he > 0 else float("nan")

        fail = ""
        if miss > a.max_site_miss:
            fail = "MISS"
        elif maf < a.min_maf:
            fail = "MAF"
        elif het > a.max_site_het:
            fail = "HET"

        if fstat:
            fstat.write(f"{chrom}\t{pos}\t{n_called}\t{miss:.4f}\t{maf:.4f}"
                        f"\t{he:.4f}\t{het:.4f}\t{fis:.4f}"
                        f"\t{0 if fail else 1}\t{fail or 'PASS'}\n")

        if fail == "MISS": drop_miss += 1; continue
        if fail == "MAF":  drop_maf  += 1; continue
        if fail == "HET":  drop_het  += 1; continue

        gts = "\t".join(CODE[c] if c in CODE else "./." for c in codes)
        fout.write(f"{chrom}\t{pos}\t{chrom}_{pos}\t{ref}\t{alt}\t.\tPASS\t"
                   f"NS={n_called};AF={q:.4f}\tGT\t{gts}\n")
        n_out += 1

    print(f"[genotype2vcf] samples kept : {len(kept_samples)}\n"
          f"[genotype2vcf] sites in     : {n_in}\n"
          f"[genotype2vcf]   dropped non-chromosome : {drop_scaf}\n"
          f"[genotype2vcf]   dropped missingness    : {drop_miss}\n"
          f"[genotype2vcf]   dropped MAF            : {drop_maf}\n"
          f"[genotype2vcf]   dropped site het       : {drop_het}\n"
          f"[genotype2vcf] sites out    : {n_out}", file=sys.stderr)

    if fout is not sys.stdout: fout.close()
    if fstat: fstat.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
