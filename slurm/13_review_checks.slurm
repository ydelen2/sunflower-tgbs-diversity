#!/bin/bash
#SBATCH --job-name=sf13_review
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=logs/sf13_review_%j.out
#SBATCH --error=logs/sf13_review_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Checks added for the first round of review: a genotype export for the
# unified diversity indices, runs of homozygosity on the non-imputed calls,
# an identity-by-state distance for permutational analysis of variance, KING
# kinship of the replicate libraries, linkage disequilibrium at a matched
# sample size, ADMIXTURE without close relatives, the replicate discordance
# breakdown and the core subset benchmark. R/24 summarizes everything.
RV=$OUT/12_review
mkdir -p "$RV" "$OUT/09_tables"
PL=$OUT/01_qc
VCF=$OUT/00_vcf/sunflower.filtered.vcf.gz           # non-imputed filtered set
VCF_IMP=$OUT/00_vcf/sunflower.filtered.imputed.vcf.gz
cd "$ANA"

# 0. list of accessions dropped for the run without close relatives
Rscript R/23_review_prepare.R

# 1. 0/1/2 genotype export of the non-imputed filtered set, all accessions
vcftools --gzvcf "$VCF" --012 --out "$RV/all"

# 2. runs of homozygosity on the non-imputed genotypes, settings of step 07
plink --bfile "$PL/sunflower" --allow-extra-chr --homozyg \
      --homozyg-density 100 --homozyg-gap 500 --homozyg-kb 500 --homozyg-snp 25 \
      --homozyg-window-snp 25 --homozyg-window-het 1 --homozyg-window-missing 5 \
      --homozyg-window-threshold 0.05 \
      --memory "$PLINK_MEM" --threads "$THREADS" --out "$RV/roh_nonimputed"

# 3. identity-by-state distance on the non-imputed set
plink --bfile "$PL/sunflower" --allow-extra-chr --distance square 1-ibs \
      --memory "$PLINK_MEM" --threads "$THREADS" --out "$RV/ibs"

# 4. KING kinship among the 276 libraries, including the three of PI 490282
plink2 --vcf "$OUT/00_vcf/mcr50.chrom.vcf.gz" --allow-extra-chr --max-alleles 2 --snps-only \
       --make-king-table --king-table-filter 0.10 \
       --memory "$PLINK_MEM" --threads "$THREADS" --out "$RV/king276"

# 5. linkage disequilibrium in ten random subsets of 23 accessions (the size of
#    Q2) drawn from the whole panel, with the settings of step 06
bcftools query -l "$VCF_IMP" > "$RV/all.ind"
for i in $(seq 1 10); do
  shuf --random-source=<(yes "$i") -n 23 "$RV/all.ind" > "$RV/sub23_$i.ind"
  PopLDdecay -InVCF "$VCF_IMP" -SubPop "$RV/sub23_$i.ind" -OutStat "$RV/sub23_$i.stat" \
             -MaxDist 1000 -MAF "$MIN_MAF" -Miss 0.2
done

# 6. ADMIXTURE at K = 5 after removing one accession of every closely related pair
DROP=$RV/drop_related.txt
if [[ -s $DROP ]]; then
  awk '{print $1"\t"$1}' "$DROP" > "$RV/drop.fam"
  plink --bfile "$PL/sunflower_pruned" --allow-extra-chr --remove "$RV/drop.fam" \
        --make-bed --memory "$PLINK_MEM" --threads "$THREADS" --out "$RV/pruned_norel"
  cd "$RV"
  for s in 1 2 3; do
    admixture -j"$THREADS" -s "$s" pruned_norel.bed 5 > "admix_norel_K5_s$s.log"
    mv pruned_norel.5.Q "pruned_norel_K5_s$s.Q"
    mv pruned_norel.5.P "pruned_norel_K5_s$s.P"
  done
  cd "$ANA"
fi

# 7. discordance between the replicate libraries by genotype class and depth
python3 scripts/replicate_discordance.py "$RAW/Sunflower.MCR50.AlleleCounts" "$PL/sunflower.bim" "$OUT/09_tables"
# read depth and allele balance at the sites removed by the heterozygous call rate ceiling
python3 scripts/site_depth_by_filter.py "$RAW/Sunflower.MCR50.AlleleCounts" "$RAW/Sunflower.MCR50.snps.genotype.txt" "$OUT/01_qc/keep_samples.txt" "$OUT/09_tables"

# 8. core subset benchmark and the summary of every check above
Rscript R/22_core_benchmark.R
Rscript R/24_review_round1.R

echo "DONE 13"
