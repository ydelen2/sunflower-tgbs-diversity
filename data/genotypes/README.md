# Filtered genotype matrices

Both files hold the 233,264 biallelic SNPs retained after the filtering of
Methods 2.3 (sites on the 17 assembled chromosomes of HanXRQr2.0-SUNRISE,
called in at least 50 percent of the retained accessions, minor allele
frequency at least 0.05, heterozygous call rate at most 0.55) for the 274
accessions of the panel. Sample names are the GRIN accession identifiers used
in Table S1. Chromosomes are named Ha1 to Ha17.

| file | content |
|---|---|
| `sunflower.filtered.vcf.gz` | genotype calls as delivered by the provider (GT only, missing calls kept); written by step 02 from `scripts/genotype2vcf.py` |
| `sunflower.filtered.imputed.vcf.gz` | the same sites and accessions taken from the provider's imputed and phased set (Beagle 5.1) |

Each file is bgzip-compressed and indexed (`.tbi`). The non-imputed file is
the input of the diversity, differentiation and nucleotide diversity steps;
the imputed file is the input of the principal component, ancestry,
discriminant, tree, kinship, runs-of-homozygosity and core subset steps and
of the main linkage disequilibrium curves (Methods 2.3).

The matrices delivered by the provider, from which these files were derived,
are available with the earlier studies of this panel at
https://github.com/ydelen2/Sunflower_important_traits (folder Genotype_data).
