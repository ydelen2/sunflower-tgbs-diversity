#!/usr/bin/env bash
# Submit the whole pipeline on Swan with the right dependencies between steps.
#   export ANA_CFG=$HOME/sunflower_tgbs/analysis/config.sh
#   bash run_all.sh
set -euo pipefail
export ANA_CFG=${ANA_CFG:-$PWD/config.sh}
source "$ANA_CFG"

cd "$ANA"

j1=$(sbatch --parsable slurm/01_make_vcf.slurm)
j2=$(sbatch --parsable --dependency=afterok:$j1 slurm/02_qc_filter.slurm)
j3=$(sbatch --parsable --dependency=afterok:$j2 slurm/03_pca.slurm)
j4=$(sbatch --parsable --dependency=afterok:$j2 slurm/04_admixture.slurm)
j7=$(sbatch --parsable --dependency=afterok:$j2 slurm/07_roh.slurm)

echo "01 make_vcf   : $j1"
echo "02 qc_filter  : $j2"
echo "03 pca        : $j3"
echo "04 admixture  : $j4  (array)"
echo "07 roh        : $j7"
cat <<'MSG'

After the ADMIXTURE array finishes:
  bash slurm/04b_admixture_collect.sh
  Rscript R/12_admixture.R          # writes results/03_admixture/group_assignment.txt
then submit the group-dependent steps:
  sbatch slurm/05_diversity_fst.slurm
  sbatch slurm/06_ld_decay.slurm
and finally the remaining R steps:
  Rscript R/10_qc_plots.R
  Rscript R/11_pca.R
  Rscript R/13_tree_dapc.R
  Rscript R/14_diversity.R
  Rscript R/15_ld_decay.R
  Rscript R/16_core_collection.R
MSG
