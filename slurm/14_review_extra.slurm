#!/bin/bash
#SBATCH --job-name=sf14_review
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sf14_review_%j.out
#SBATCH --error=logs/sf14_review_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Added after slurm/13: the unfiltered chromosome-anchored calls of the retained
# accessions, so that the sites removed by the heterozygous call rate filter can
# be examined within each cluster (R/24, section N). Then R/24 is run again so
# that every review table is current.
RV=$OUT/12_review
mkdir -p "$RV"
cd "$ANA"
vcftools --gzvcf "$OUT/00_vcf/mcr50.chrom.vcf.gz" --keep "$OUT/01_qc/keep_samples.txt" --012 --out "$RV/allsites"
Rscript R/24_review_round1.R
# the linkage disequilibrium figure gains a ratio panel that reads the subset curves written by R/24
Rscript R/20_redraw_figures.R
echo "DONE 14"
