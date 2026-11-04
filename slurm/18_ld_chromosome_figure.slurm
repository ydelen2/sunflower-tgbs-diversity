#!/bin/bash
#SBATCH --job-name=sf18_ld_chr_fig
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/sf18_ld_chr_fig_%j.out
#SBATCH --error=logs/sf18_ld_chr_fig_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Figure S10: the whole-panel linkage disequilibrium decay of every chromosome
# (from the per-chromosome runs of slurm/12) and the value at 100 kb per chromosome.
cd "$ANA"
Rscript R/26_ld_by_chromosome_figure.R
ls -la "$OUT/08_figures/FigS10_ld_by_chromosome."*
echo "DONE 18"
