#!/bin/bash
#SBATCH --job-name=sf16_masksum
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sf16_masksum_%j.out
#SBATCH --error=logs/sf16_masksum_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Summary step of the masking experiment (slurm/15), rerun on its own: the
# genotype tables are rewritten as integer dosages and R/25 is run again.
MK=$OUT/12_review/masking
cd "$ANA"
for tag in mask0 mask10 mask20; do
  bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$MK/$tag.imputed.vcf.gz" \
    | sed -e 's#0|0#0#g' -e 's#1|1#2#g' -e 's#0|1#1#g' -e 's#1|0#1#g' | gzip -c > "$MK/$tag.imputed.gt.tsv.gz"
done
Rscript R/25_masking_summary.R
echo "DONE 16"
