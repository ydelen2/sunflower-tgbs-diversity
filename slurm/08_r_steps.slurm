#!/bin/bash
#SBATCH --job-name=sf08_R
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sf08_R_%j.out
#SBATCH --error=logs/sf08_R_%j.err
set -euo pipefail

source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Runs the R analysis steps. The tree bootstrap and the AMOVA permutation test
# are too heavy for a login node at 274 accessions and 49,523 pruned markers.
# Pass a subset as arguments to rerun only those, for example:
#   sbatch slurm/08_r_steps.slurm 13_tree_dapc 14_diversity

cd "$ANA"
STEPS=("$@")
if [[ ${#STEPS[@]} -eq 0 ]]; then
  STEPS=(10_qc_plots 11_pca 12_admixture 13_tree_dapc 14_diversity 15_ld_decay 16_core_collection 17_manuscript_stats 18_fst_intervals 19_supplementary_additions 20_redraw_figures)
fi

for s in "${STEPS[@]}"; do
  echo "===== $s ====="
  Rscript "R/${s}.R"
  echo "===== $s finished ====="
done

echo "DONE 08"
