#!/bin/bash
#SBATCH --job-name=sf05_div
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sf05_div_%j.out
#SBATCH --error=logs/sf05_div_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

DV=$OUT/04_diversity
VCF=$OUT/00_vcf/sunflower.filtered.vcf.gz   # non-imputed: diversity statistics are not
                                            # estimated from Beagle-filled genotypes
GRP=$OUT/03_admixture/group_assignment_pops.txt   # written by R/12_admixture.R:
                                            # SAMPLE, GROUP; clusters only

cd "$DV"

# ---------------------------------------------------------------------------
# 1. Whole-collection statistics
# ---------------------------------------------------------------------------
vcftools --gzvcf "$VCF" --het            --out all      # per-accession F
vcftools --gzvcf "$VCF" --site-pi        --out all
vcftools --gzvcf "$VCF" --window-pi 1000000 --window-pi-step 200000 --out all
vcftools --gzvcf "$VCF" --TajimaD 1000000               --out all
vcftools --gzvcf "$VCF" --freq2          --out all
vcftools --gzvcf "$VCF" --hardy          --out all
vcftools --gzvcf "$VCF" --SNPdensity 1000000 --out all

# ---------------------------------------------------------------------------
# 2. Per-group statistics + pairwise Weir & Cockerham Fst
# ---------------------------------------------------------------------------
if [[ -s $GRP ]]; then
  mkdir -p pops
  cut -f2 "$GRP" | sort -u > pops/group_names.txt
  while read -r g; do
      awk -v G="$g" -F'\t' '$2==G {print $1}' "$GRP" > "pops/${g}.ind"
  done < pops/group_names.txt

  while read -r g; do
      n=$(wc -l < "pops/${g}.ind")
      echo "== group $g (n=$n) =="
      [[ $n -lt 4 ]] && { echo "   skipped: n<4"; continue; }
      vcftools --gzvcf "$VCF" --keep "pops/${g}.ind" --het               --out "grp_${g}"
      vcftools --gzvcf "$VCF" --keep "pops/${g}.ind" --site-pi           --out "grp_${g}"
      vcftools --gzvcf "$VCF" --keep "pops/${g}.ind" --window-pi 1000000 --window-pi-step 200000 --out "grp_${g}"
      vcftools --gzvcf "$VCF" --keep "pops/${g}.ind" --TajimaD 1000000   --out "grp_${g}"
      vcftools --gzvcf "$VCF" --keep "pops/${g}.ind" --freq2             --out "grp_${g}"
  done < pops/group_names.txt

  # pairwise Fst, all group pairs
  mapfile -t G < pops/group_names.txt
  for ((i=0; i<${#G[@]}; i++)); do
    for ((j=i+1; j<${#G[@]}; j++)); do
      a=${G[$i]}; b=${G[$j]}
      [[ $(wc -l < "pops/${a}.ind") -lt 4 || $(wc -l < "pops/${b}.ind") -lt 4 ]] && continue
      # vcftools writes the genome-wide Fst estimates to stderr, not to a .log
      # file, so each call's stderr is captured per pair and parsed below.
      vcftools --gzvcf "$VCF" \
               --weir-fst-pop "pops/${a}.ind" --weir-fst-pop "pops/${b}.ind" \
               --fst-window-size 100000 --fst-window-step 100000 \
               --out "fst_${a}_vs_${b}" 2> "fst_${a}_vs_${b}.vcftools.log"
      # per-site Fst values
      vcftools --gzvcf "$VCF" \
               --weir-fst-pop "pops/${a}.ind" --weir-fst-pop "pops/${b}.ind" \
               --out "fstsite_${a}_vs_${b}" 2> "fstsite_${a}_vs_${b}.vcftools.log"
    done
  done

  { printf 'POP1\tPOP2\tN1\tN2\tMEAN_FST\tWEIGHTED_FST\n'
    for f in fst_*.vcftools.log; do
      [[ -e $f ]] || continue
      p=$(basename "$f" .vcftools.log); p=${p#fst_}
      a=${p%%_vs_*}; b=${p##*_vs_}
      m=$(awk '/Weir and Cockerham mean Fst estimate/{v=$NF} END{print v}'     "$f")
      w=$(awk '/Weir and Cockerham weighted Fst estimate/{v=$NF} END{print v}' "$f")
      printf '%s\t%s\t%d\t%d\t%s\t%s\n' "$a" "$b" \
             "$(wc -l < "pops/${a}.ind")" "$(wc -l < "pops/${b}.ind")" \
             "${m:-NA}" "${w:-NA}"
    done; } > pairwise_fst_summary.tsv

  cat pairwise_fst_summary.tsv
else
  echo "NOTE: $GRP not found - run R/12_admixture.R first, then re-run this step."
fi

echo "DONE 05"
