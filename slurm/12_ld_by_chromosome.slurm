#!/bin/bash
#SBATCH --job-name=sf12_ldchr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/sf12_ldchr_%j.out
#SBATCH --error=logs/sf12_ldchr_%j.err
set -euo pipefail

source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Whole-panel LD decay computed for each chromosome separately, with the same
# settings as slurm/06_ld_decay.slurm (imputed filtered set, pairs within 1 Mb,
# within-set MAF 0.05). Output feeds R/19_supplementary_additions.R.

LD=$OUT/05_ld/bychr
VCF=$OUT/00_vcf/sunflower.filtered.imputed.vcf.gz
mkdir -p "$LD"
cd "$LD"

[[ -s $VCF.tbi || -s $VCF.csi ]] || bcftools index -f "$VCF"

MAXDIST=1000        # kb, as in slurm/06
while read -r NC N; do
  [[ -n "${NC:-}" ]] || continue
  # the filtered VCF names chromosomes 1 to 17 (second column of CHRMAP)
  bcftools view -r "$N" -Oz -o "chr_Ha${N}.vcf.gz" "$VCF"
  PopLDdecay -InVCF "chr_Ha${N}.vcf.gz" -OutStat "chr_Ha${N}.stat" \
             -MaxDist $MAXDIST -MAF "$MIN_MAF" -Miss 0.2
  rm -f "chr_Ha${N}.vcf.gz"
  echo "chromosome Ha${N} done"
done < "$CHRMAP"

ls -1 chr_Ha*.stat* | wc -l
echo "DONE 12"
