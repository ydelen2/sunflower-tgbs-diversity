#!/usr/bin/env python3
"""Mask a random fraction of the called genotypes in a VCF, for a masking
experiment on imputation. Every called genotype of every sample is set to
missing with the given probability; the masked cells and their original
genotypes are written to a tab-separated file so that the imputed values can
be scored afterwards.

usage: mask_genotypes.py IN_VCF_GZ FRACTION SEED OUT_VCF_GZ MASKED_TSV_GZ
"""
import gzip
import random
import sys


def main(vcf_in, fraction, seed, vcf_out, masked_out):
    fraction = float(fraction)
    rng = random.Random(int(seed))
    n_called = n_masked = 0
    samples = []
    with gzip.open(vcf_in, "rt") as fin, gzip.open(vcf_out, "wt") as fout, gzip.open(masked_out, "wt") as fm:
        fm.write("CHROM\tPOS\tSAMPLE\tGT\n")
        for line in fin:
            if line.startswith("##"):
                fout.write(line); continue
            if line.startswith("#"):
                fout.write(line); samples = line.rstrip("\n").split("\t")[9:]; continue
            f = line.rstrip("\n").split("\t")
            fmt = f[8].split(":")
            gi = fmt.index("GT")
            for k in range(9, len(f)):
                parts = f[k].split(":")
                gt = parts[gi]
                if gt in ("./.", ".", ".|."):
                    continue
                n_called += 1
                if fraction > 0 and rng.random() < fraction:
                    fm.write("%s\t%s\t%s\t%s\n" % (f[0], f[1], samples[k - 9], gt.replace("|", "/")))
                    parts[gi] = "./."
                    f[k] = ":".join(parts)
                    n_masked += 1
            fout.write("\t".join(f) + "\n")
    print("called genotypes %d, masked %d (%.4f)" % (n_called, n_masked, n_masked / max(n_called, 1)))
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:6]))
