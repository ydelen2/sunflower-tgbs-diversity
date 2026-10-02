#!/bin/bash
#SBATCH --job-name=sf07_roh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sf07_roh_%j.out
#SBATCH --error=logs/sf07_roh_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Runs of homozygosity -> genomic inbreeding coefficient F_ROH.
# At the marker density of this set (233,264 SNPs, 1 per 12.9 kb) a 500-kb
# segment holds about 39 SNPs, so --homozyg-kb 500 is the binding criterion
# and --homozyg-snp 25 is not. Minimum reported ROH length is 500 kb.
plink --bfile "$OUT/01_qc/sunflower_imp" --allow-extra-chr \
      --homozyg \
      --homozyg-density 100 \
      --homozyg-gap 500 \
      --homozyg-kb 500 \
      --homozyg-snp 25 \
      --homozyg-window-snp 25 \
      --homozyg-window-het 1 \
      --homozyg-window-missing 5 \
      --homozyg-window-threshold 0.05 \
      --memory "$PLINK_MEM" --threads "$THREADS" \
      --out "$OUT/04_diversity/roh"

echo "DONE 07"
