#!/bin/bash
# Collect ADMIXTURE cross-validation errors and log-likelihoods across the grid.
set -euo pipefail
source "${ANA_CFG:-$HOME/sunflower_tgbs/analysis/config.sh}"

AD=$OUT/03_admixture
{
  echo -e "K\tREP\tCV_ERROR\tLOGLIK"
  for d in "$AD"/K*_rep*; do
      [[ -d $d ]] || continue
      k=$(basename "$d" | sed 's/K\([0-9]*\)_rep.*/\1/')
      r=$(basename "$d" | sed 's/.*_rep//')
      lg=$(ls "$d"/log_K*_rep*.out 2>/dev/null | head -1) || continue
      cv=$(grep -o 'CV error (K=[0-9]*): [0-9.]*' "$lg" | awk '{print $NF}')
      ll=$(grep -E '^Loglikelihood:' "$lg" | tail -1 | awk '{print $2}')
      echo -e "${k}\t${r}\t${cv:-NA}\t${ll:-NA}"
  done
} > "$AD/cv_summary.tsv"

# best (lowest CV) replicate per K -> used for the barplots
awk 'NR>1 && $3!="NA"' "$AD/cv_summary.tsv" | sort -k1,1n -k3,3g \
  | awk '!seen[$1]++ {print $1"\t"$2"\t"$3}' > "$AD/best_rep_per_K.tsv"

cat "$AD/best_rep_per_K.tsv"
echo "wrote $AD/cv_summary.tsv"
