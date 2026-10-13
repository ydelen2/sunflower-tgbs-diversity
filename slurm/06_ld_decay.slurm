#!/bin/bash
#SBATCH --job-name=sf06_ld
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sf06_ld_%j.out
#SBATCH --error=logs/sf06_ld_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

LD=$OUT/05_ld
VCF=$OUT/00_vcf/sunflower.filtered.imputed.vcf.gz   # imputed set: every accession
                                                    # contributes to every pair. The non-imputed
                                                    # curve is compared in slurm/11_ld_sensitivity
                                                    # and Table S16.
GRP=$OUT/03_admixture/group_assignment_pops.txt
cd "$LD"

MAXDIST=1000        # kb - sunflower LD in cultivated material can extend far
BIN=1               # kb

# ---------------------------------------------------------------------------
# Option A (preferred): PopLDdecay
# ---------------------------------------------------------------------------
if command -v PopLDdecay >/dev/null 2>&1; then
  PopLDdecay -InVCF "$VCF" -OutStat all.stat -MaxDist $MAXDIST -MAF "$MIN_MAF" -Miss 0.2
  if [[ -s $GRP ]]; then
    cut -f2 "$GRP" | sort -u | while read -r g; do
      awk -v G="$g" -F'\t' '$2==G {print $1}' "$GRP" > "sub_${g}.ind"
      [[ $(wc -l < "sub_${g}.ind") -ge 20 ]] || { echo "skip $g (n<20)"; continue; }
      PopLDdecay -InVCF "$VCF" -SubPop "sub_${g}.ind" -OutStat "grp_${g}.stat" \
                 -MaxDist $MAXDIST -MAF "$MIN_MAF" -Miss 0.2
    done
  fi

# ---------------------------------------------------------------------------
# Option B (fallback, no extra software): plink pairwise r2 + binning
# ---------------------------------------------------------------------------
else
  echo "PopLDdecay not found - falling back to plink --r2"
  for CHR in $(cut -f2 "$CHRMAP"); do
    plink --bfile "$OUT/01_qc/sunflower_imp" --allow-extra-chr \
          --chr "$CHR" --r2 gz --ld-window 999999 \
          --ld-window-kb $MAXDIST --ld-window-r2 0 \
          --memory "$PLINK_MEM" --threads "$THREADS" \
          --out "chr${CHR}"
  done
  # bin into 1-kb classes: distance -> mean r2, n pairs
  zcat chr*.ld.gz \
    | awk -v b=$((BIN*1000)) 'NR>1 && $1!="CHR_A"{d=$5-$2; if(d<0)d=-d; k=int(d/b);
        s[k]+=$7; c[k]++} END{ printf "DIST_KB\tMEAN_R2\tN_PAIRS\n";
        for(k in s) printf "%.1f\t%.6f\t%d\n", (k+0.5)*b/1000, s[k]/c[k], c[k]}' \
    | sort -k1,1g > all.bin.tsv
fi

echo "DONE 06"
