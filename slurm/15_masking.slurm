#!/bin/bash
#SBATCH --job-name=sf15_mask
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sf15_mask_%j.out
#SBATCH --error=logs/sf15_mask_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Masking experiment on imputation and runs of homozygosity (peer review,
# comment 16). The non-imputed filtered calls of the 274 accessions are
# re-imputed with Beagle three times: as delivered (no masking), and after
# masking a random 10 and 20 percent of the called genotypes. The imputed
# genotypes at the masked cells are scored against the originals, and runs of
# homozygosity are called on each imputed set with the settings of step 07.
# R/25 summarizes the accuracy and compares F_ROH with the provider's imputation.
MK=$OUT/12_review/masking
mkdir -p "$MK" "$PROJ/tools"
cd "$ANA"
BEAGLE=$PROJ/tools/beagle.22Jul22.46e.jar
if [[ ! -s "$BEAGLE" ]]; then
  wget -q -O "$BEAGLE" https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
fi
java -version 2>&1 | head -1
VCF=$OUT/00_vcf/sunflower.filtered.vcf.gz

for frac in 0 0.10 0.20; do
  tag=mask${frac/0./}
  python3 scripts/mask_genotypes.py "$VCF" "$frac" 2026 "$MK/$tag.vcf.gz" "$MK/$tag.masked.tsv.gz"
  java -Xmx48g -jar "$BEAGLE" gt="$MK/$tag.vcf.gz" out="$MK/$tag.imputed" nthreads="$THREADS" seed=2026
  # genotypes as integer dosages (0, 1, 2) so that the summary step stays small
  bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$MK/$tag.imputed.vcf.gz" \
    | sed -e 's#0|0#0#g' -e 's#1|1#2#g' -e 's#0|1#1#g' -e 's#1|0#1#g' | gzip -c > "$MK/$tag.imputed.gt.tsv.gz"
  plink --vcf "$MK/$tag.imputed.vcf.gz" --allow-extra-chr --double-id \
        --memory "$PLINK_MEM" --threads "$THREADS" --make-bed --out "$MK/$tag"
  plink --bfile "$MK/$tag" --allow-extra-chr \
        --homozyg --homozyg-density 100 --homozyg-gap 500 --homozyg-kb 500 --homozyg-snp 25 \
        --homozyg-window-snp 25 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 \
        --memory "$PLINK_MEM" --threads "$THREADS" --out "$MK/roh_$tag"
done
bcftools query -l "$VCF" > "$MK/samples.txt"
Rscript R/25_masking_summary.R
echo "DONE 15"
