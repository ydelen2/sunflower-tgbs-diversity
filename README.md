# Sunflower tGBS diversity and population structure

Analysis code and derived data for a panel of 274 cultivated sunflower
(*Helianthus annuus* L.) accessions from the USDA National Plant Germplasm
System, genotyped with tunable genotyping by sequencing and analyzed against
the HanXRQr2.0-SUNRISE assembly (GenBank GCF_002127325.2).

## Repository layout

| path | content |
|---|---|
| `config.sh` | paths, filtering thresholds and the ADMIXTURE grid; every script sources it |
| `run_all.sh` | submits the cluster steps in dependency order |
| `slurm/` | batch steps 01 to 17 |
| `R/` | analysis and figure steps 00 to 25 |
| `scripts/` | genotype conversion, per-accession statistics, GRIN retrieval, multi-panel figure assembly, replicate library discordance, read depth by site class, genotype masking, kinship-corrected linkage disequilibrium |
| `env/environment.yml` | conda environment |
| `ref/` | chromosome names and lengths of the reference assembly |
| `metadata/` | GRIN passport, descriptor and type records, sample lists, per-accession quality summary |
| `data/` | nine large per-accession or per-window result files, and under `data/genotypes/` the filtered genotype matrices (VCF, non-imputed and imputed) |
| `tables/main/` | the three main tables of the manuscript |
| `tables/supplementary/` | the 16 supplementary tables; Tables S3, S5, S6, S7, S10, S11, S13, S14, S15 and S16 hold further blocks, in the files ending in `b` to `g` |
| `tables/additional/` | replicate concordance, duplicate pairs and the imputation check of the principal components, whose values are given in the text |

## Input not included here

The genotype delivery (247,008 and 503,188 site matrices, the imputed VCF and
the per-accession allele counts) and the filtered VCF are too large for this
repository and are archived separately; the accession number is given in the
data availability statement of the associated manuscript. Genotyping used
Bsp1286I, an Illumina HiSeq X instrument and GSNAP alignment, with only
uniquely aligned reads used for variant calling.

## Filtering thresholds

| parameter | value | applies to |
|---|---|---|
| `MAX_SAMPLE_MISS` | 0.90 | accession missing rate |
| `MAX_SITE_MISS` | 0.50 | site call rate |
| `MIN_MAF` | 0.05 | minor allele frequency |
| `MAX_SITE_HET` | 0.55 | heterozygous call rate, removes collapsed paralogs |
| `LD_PRUNE` | 50 10 0.2 | `plink --indep-pairwise` window, step, r2 |
| `KING_DUP` | 0.354 | kinship above which a pair is treated as one accession |
| `KMIN` `KMAX` `NREP` `ADMIX_CV` | 1 12 10 10 | ADMIXTURE grid |

Filtering retains 274 accessions and 233,264 SNPs, of which 49,523 remain after
pruning for linkage disequilibrium.

## Order of the cluster steps

```bash
export PROJ=/path/to/project
export ANA_CFG=$PROJ/analysis/config.sh
cd $PROJ/analysis && bash run_all.sh
```

| step | script | main output |
|---|---|---|
| 01 | `slurm/01_make_vcf.slurm` | chromosome-anchored VCF, imputed VCF, per-site and per-accession statistics |
| 02 | `slurm/02_qc_filter.slurm` | filtered VCF, plink binaries, kinship, genomic relationship matrix, pruned set |
| 03 | `slurm/03_pca.slurm` | principal components on imputed and non-imputed genotypes |
| 04 | `slurm/04_admixture.slurm` | ADMIXTURE runs, K = 1 to 12, ten seeds per K |
| 04b | `slurm/04b_admixture_collect.sh` | cross-validation summary, best run per K |
| 05 | `slurm/05_diversity_fst.slurm` | per-cluster diversity, nucleotide diversity, pairwise F_ST |
| 06 | `slurm/06_ld_decay.slurm` | linkage disequilibrium decay, whole panel and per cluster |
| 07 | `slurm/07_roh.slurm` | runs of homozygosity |
| 08 | `slurm/08_r_steps.slurm` | runs the R steps that exceed login node limits |
| 09 | `slurm/09_allsnp_comparison.slurm` | call rate threshold comparison on the 503,188 site delivery |
| 10 | `slurm/10_region_fst.slurm` | pairwise F_ST among geographic regions |
| 11 | `slurm/11_ld_sensitivity.slurm` | linkage disequilibrium without imputation and under a stricter heterozygosity filter |
| 12 | `slurm/12_ld_by_chromosome.slurm` | linkage disequilibrium decay per chromosome |
| 13 | `slurm/13_review_checks.slurm` | checks added at peer review: relatedness pruning, non-imputed runs of homozygosity, identity-by-state distances, KING on the 276 libraries, linkage disequilibrium in ten random subsets of 23 accessions, ADMIXTURE without close relatives, replicate discordance, core benchmark; then `R/22` to `R/24` |
| 14 | `slurm/14_review_extra.slurm` | unfiltered calls of the retained accessions for the removed-site check; reruns `R/24` and `R/20` |
| 15 | `slurm/15_masking.slurm` | masking experiment: the non-imputed calls re-imputed with Beagle 5.4 as delivered and after masking 10 and 20 percent of the called genotypes, with runs of homozygosity on each set; then `R/25` |
| 16 | `slurm/16_masking_summary.slurm` | the summary step of 15 on its own |
| 17 | `slurm/17_depth_and_ld.slurm` | read depth of the removed and retained sites (Table S3f) and linkage disequilibrium corrected for kinship (Table S14d), both from the delivered matrices |

## Order of the R steps

| script | reads | writes |
|---|---|---|
| `R/01_build_metadata.R` | GRIN records | `metadata/sample_metadata.csv` |
| `R/10_qc_plots.R` | steps 01, 02 | Table S4, Table S3, Data file 2, replicate concordance |
| `R/11_pca.R` | steps 02, 03 | Table S8, Data file 4 |
| `R/12_admixture.R` | step 04b | Tables S5, S6 and S7, Data file 6, cluster assignment |
| `R/13_tree_dapc.R` | step 02, `R/12` | Tables S2 and S7, Data file 3 |
| `R/14_diversity.R` | steps 05, 07 | Table 1, Tables S10, S11 and S12, Data files 7 and 8, pairwise F_ST used by `R/18` |
| `R/15_ld_decay.R` | steps 06, 11 | Table 2, Table S16, Data file 9 |
| `R/16_core_collection.R` | step 02, `R/12` | Table 3, Table S15, Data file 1 |
| `R/17_manuscript_stats.R` | steps 03, `R/11` to `R/14` | Tables S1, S5, S9 and S13, imputation check of the principal components |
| `R/18_fst_intervals.R` | step 02, `R/12`, `R/14` | Table S10 |
| `R/19_supplementary_additions.R` | step 12, `R/16`, `R/17` | Tables S14 and S15 |
| `R/20_redraw_figures.R` | results of `R/10` to `R/16`, `R/24` | every figure except Figure S4, as published, including panel (c) of Figure 4 |
| `R/21_core_repeat.R` | step 02, `R/16` | Table S15 (repeat of the core subset searches) |
| `R/22_core_benchmark.R` | step 02, `R/16` | Table S15d (the two cores against 100 random and 100 proportional stratified random subsets) |
| `R/23_review_prepare.R` | step 02 | Table S6b (relatedness within and between clusters) and the list of close relatives removed in step 13 |
| `R/24_review_round1.R` | steps 02 to 07, 13 and 14, `R/12`, `R/16` | the tables added at peer review: Table 1 (unified diversity indices, receipt years), Table 2 (random subsets), Tables S1, S9, S10b, S15b (extended), S3b to S3e, S5c, S6c, S7c, S10d to S10f, S11d, S13b to S13e, S14c, S16d |

| `R/25_masking_summary.R` | step 15, step 07 | Tables S13f and S13g (accuracy at the masked cells, genomic inbreeding of the re-imputed sets) |

`scripts/ld_kinship_corrected.py` computes linkage disequilibrium corrected for
kinship and for kinship and ancestry (Mangin et al. 2012) on the imputed
genotype matrix, with the genomic relationship matrix (`data/Data1_kinship.csv`),
Tables S1 and S6 and the list of retained sites as inputs, and writes
Table S14d; step 17 runs it.
`scripts/replicate_discordance.py` compares the three libraries of the control
accession on the delivered allele counts and writes Tables S3b and S3c, and
`scripts/site_depth_by_filter.py` compares read depth and allele balance at the
sites removed by the heterozygous call rate ceiling with the retained sites
(Table S3f); step 13 calls both, and step 17 repeats the depth comparison before the kinship-corrected linkage disequilibrium.

`R/12_admixture.R` writes the cluster assignment of every accession, including
those below the ancestry threshold of 0.60, which are labeled admixed. Steps
05 and 06 read the file that holds the assigned clusters only.

## Data files

| file | content |
|---|---|
| `data/Data1_kinship.csv` | genomic relationship matrix of the 274 accessions |
| `data/Data2_accession_qc.csv` | per-accession call rate, missing rate and heterozygosity |
| `data/Data3_dapc_coordinates.csv` | discriminant function coordinates, group and cluster |
| `data/Data4_pca_scores.csv` | scores on the first twenty principal components with passport data |
| `data/Data5_bionj_tree_rect.pdf` | the distance tree with tip labels and bootstrap support |
| `data/Data6_K3_assignment.csv` | cluster assignment at K = 3 with the maximum ancestry coefficient |
| `data/Data7_ROH.csv` | runs of homozygosity per accession |
| `data/Data8_windowed_pi.csv` | nucleotide diversity per 1-Mb window and cluster |
| `data/Data9_LD_binned_empirical.csv` | mean r-squared and marker pair counts in 2-kb bins |
| `data/genotypes/sunflower.filtered.vcf.gz` | the filtered genotype calls, 233,264 sites by 274 accessions, non-imputed, with tabix index |
| `data/genotypes/sunflower.filtered.imputed.vcf.gz` | the same sites and accessions from the imputed and phased set of the provider (Beagle 5.1), with tabix index |

## Software

`env/environment.yml` holds the conda environment. Two tools are installed
separately:

```bash
R -e 'install.packages("corehunter", repos="https://cloud.r-project.org")'
git clone https://github.com/BGI-shenzhen/PopLDdecay && cd PopLDdecay && ./configure && make
```

Versions used: PLINK 1.9.0-b.8 and PLINK 2.0.0-a.6.9LM, VCFtools 0.1.17,
BCFtools 1.24 with HTSlib 1.24, ADMIXTURE 1.3.0, PopLDdecay commit 5dcd070,
Beagle 5.1 (applied by the genotyping provider), R 4.3.3 with SNPRelate 1.36.0,
adegenet 2.1.11, ape 5.8.1, poppr 2.9.8, hierfstat 0.5.11, vcfR 1.15.0,
corehunter 3.2.3, data.table 1.17.8, ggtree 3.14.0 and ggplot2 3.5.2, and
Python 3.11 with numpy 2.4 and pandas 3.1.

## License

Code is released under the MIT license, see `LICENSE`. The files under `data/`,
`tables/` and `metadata/` are released under CC BY 4.0,
https://creativecommons.org/licenses/by/4.0/. Passport and descriptor records
originate from the GRIN-Global database of the USDA National Plant Germplasm
System.
