#!/bin/bash
#SBATCH --job-name=sf02_qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/sf02_qc_%j.out
#SBATCH --error=logs/sf02_qc_%j.err
set -euo pipefail

source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

QC=$OUT/01_qc
VC=$OUT/00_vcf
cd "$VC"

# ---------------------------------------------------------------------------
# 1. Accession keep-list: missingness threshold, minus the genotyping-control
#    replicates.  PI_490282 was added once per plate as a plate control
#    (PI_490282C2, PI_490282C3); they enter the reproducibility check below but
#    must not enter population analyses.
# ---------------------------------------------------------------------------
awk -v t="$MAX_SAMPLE_MISS" 'NR>1 && $4<=t {print $1}' "$QC/sample_stats.tsv" \
  | grep -v -E '^PI_490282C[0-9]+$' > "$QC/keep_samples.txt"
awk -v t="$MAX_SAMPLE_MISS" 'NR>1 && $4> t {print $1"\t"$4}' "$QC/sample_stats.tsv" > "$QC/dropped_samples.txt"

# genotyping reproducibility across the three PI_490282 libraries, computed on
# the full MCR50 matrix
awk -F'\t' '
NR==2{for(i=5;i<=NF;i++){if($i=="PI_490282")a=i; if($i=="PI_490282C2")b=i; if($i=="PI_490282C3")c=i} next}
NR>2{ if($a!=0&&$b!=0){t1++; if($a==$b)m1++}
      if($a!=0&&$c!=0){t2++; if($a==$c)m2++}
      if($b!=0&&$c!=0){t3++; if($b==$c)m3++} }
END{printf "PAIR\tCONCORDANCE\tN_BOTH_CALLED\nPI_490282_vs_C2\t%.4f\t%d\nPI_490282_vs_C3\t%.4f\t%d\nC2_vs_C3\t%.4f\t%d\n", m1/t1,t1, m2/t2,t2, m3/t3,t3}' \
  "$GENO_MCR50" > "$QC/control_replicate_concordance.tsv"
cat "$QC/control_replicate_concordance.tsv"
N_COL=$(( $(wc -l < "$QC/sample_stats.tsv") - 1 ))
N_CTRL=$(awk 'NR>1 && $1 ~ /^PI_490282C[0-9]+$/' "$QC/sample_stats.tsv" | wc -l)
N_DROP=$(wc -l < "$QC/dropped_samples.txt")
N_KEEP=$(wc -l < "$QC/keep_samples.txt")
{
  printf 'CATEGORY\tN\n'
  printf 'columns in the MCR50 matrix\t%d\n' "$N_COL"
  printf 'plate-control replicate libraries\t%d\n' "$N_CTRL"
  printf 'dropped, missingness above %s\t%d\n' "$MAX_SAMPLE_MISS" "$N_DROP"
  printf 'retained for analysis\t%d\n' "$N_KEEP"
} > "$QC/accession_accounting.tsv"
cat "$QC/accession_accounting.tsv"

# ---------------------------------------------------------------------------
# 2. Rebuild the working VCF on the retained accessions, applying site filters.
#    Order matters: site call rate must be recomputed AFTER removing bad
#    accessions, otherwise a handful of poor libraries drives the whole set.
# ---------------------------------------------------------------------------
python3 "$ANA/scripts/genotype2vcf.py" \
    -i "$GENO_MCR50" \
    --chr-map "$CHRMAP" --chrom-only \
    --keep "$QC/keep_samples.txt" \
    --max-site-miss "$MAX_SITE_MISS" \
    --min-maf "$MIN_MAF" \
    --max-site-het "$MAX_SITE_HET" \
    --stats "$QC/site_stats_filtered.tsv" \
  | bgzip -@ "$THREADS" -c > sunflower.filtered.vcf.gz
tabix -f -p vcf sunflower.filtered.vcf.gz
# matching subset of the imputed set (same accessions, same sites)
bcftools view -S "$QC/keep_samples.txt" --force-samples \
    -R <(bcftools query -f '%CHROM\t%POS\n' sunflower.filtered.vcf.gz) \
    mcr50.imputed.chrom.vcf.gz -Oz -o sunflower.filtered.imputed.vcf.gz --threads "$THREADS"
tabix -f -p vcf sunflower.filtered.imputed.vcf.gz

# ---------------------------------------------------------------------------
# 2b. Recompute the accession QC table on the retained set only.  Site call
#     rates, and therefore the "common well-covered sites" used for the
#     depth-corrected Ho, improve substantially once the poor libraries are out.
#     This table supplies the CORE_HO column reported in the Methods.
# ---------------------------------------------------------------------------
python3 "$ANA/scripts/sample_stats.py" \
    -i "$GENO_MCR50" --chr-map "$CHRMAP" --core-call 0.90 \
    --keep "$QC/keep_samples.txt" \
    -o "$QC/sample_stats_kept.tsv"

# ---------------------------------------------------------------------------
# 3. PLINK binaries
# ---------------------------------------------------------------------------
PL=$OUT/01_qc
plink2 --vcf sunflower.filtered.vcf.gz \
       --double-id --allow-extra-chr --set-missing-var-ids '@_#' \
       --max-alleles 2 --snps-only \
       --memory "$PLINK_MEM" --threads "$THREADS" \
       --make-bed --out "$PL/sunflower"

plink2 --vcf sunflower.filtered.imputed.vcf.gz \
       --double-id --allow-extra-chr --set-missing-var-ids '@_#' \
       --max-alleles 2 --snps-only \
       --memory "$PLINK_MEM" --threads "$THREADS" \
       --make-bed --out "$PL/sunflower_imp"

# ---------------------------------------------------------------------------
# 4. Duplicate / clone detection (KING-robust kinship on the imputed set)
#    0.354 ~ identical genotype; 0.177 ~ full sib / parent-offspring
# ---------------------------------------------------------------------------
plink2 --bfile "$PL/sunflower_imp" --allow-extra-chr \
       --make-king-table --king-table-filter 0.15 \
       --memory "$PLINK_MEM" --threads "$THREADS" \
       --out "$PL/king_pairs"
awk -v t="$KING_DUP" 'NR==1 || $NF>=t' "$PL/king_pairs.kin0" > "$PL/king_duplicates.tsv"

# Candidate duplicates are judged against the replicate baseline by direct
# genotype concordance, not by the KING coefficient alone: KING assumes outbred
# samples and is inflated in inbred material.
{
  printf 'ACCESSION1\tACCESSION2\tKING\tCONCORDANCE\tN_BOTH_CALLED\n'
  awk 'NR>1 {print $2, $4, $NF}' "$PL/king_duplicates.tsv" | while read -r a b k; do
    awk -F'\t' -v A="$a" -v B="$b" -v K="$k" '
      NR==2 { for (i=5;i<=NF;i++) { if ($i==A) x=i; if ($i==B) y=i } next }
      NR>2 && x && y && $x!=0 && $y!=0 { n++; if ($x==$y) m++ }
      END { if (n) printf "%s\t%s\t%s\t%.4f\t%d\n", A, B, K, m/n, n }' "$GENO_MCR50"
  done
} > "$QC/duplicate_pairs.tsv"
cat "$QC/duplicate_pairs.tsv"
mkdir -p "$OUT/09_tables"
cp "$QC/duplicate_pairs.tsv" "$OUT/09_tables/duplicate_pairs.tsv"

# genomic relationship matrix, GCTA form (Yang et al. 2011), for the kinship heat map
plink2 --bfile "$PL/sunflower_imp" --allow-extra-chr \
       --make-grm-list --memory "$PLINK_MEM" --threads "$THREADS" \
       --out "$PL/grm"

# ---------------------------------------------------------------------------
# 5. LD-pruned marker set for PCA / ADMIXTURE / NJ tree
# ---------------------------------------------------------------------------
# shellcheck disable=SC2086
plink2 --bfile "$PL/sunflower_imp" --allow-extra-chr \
       --indep-pairwise $LD_PRUNE \
       --memory "$PLINK_MEM" --threads "$THREADS" \
       --out "$PL/prune"

plink2 --bfile "$PL/sunflower_imp" --allow-extra-chr \
       --extract "$PL/prune.prune.in" \
       --memory "$PLINK_MEM" --threads "$THREADS" \
       --make-bed --out "$PL/sunflower_pruned"

# additive 0/1/2 export, consumed by R/16_core_collection.R
plink2 --bfile "$PL/sunflower_pruned" --allow-extra-chr \
       --export A --memory "$PLINK_MEM" --threads "$THREADS" \
       --out "$PL/sunflower_pruned_A"

# ADMIXTURE refuses non-integer chromosome names; ours are already 1..17
awk '{$1=$1; print}' OFS='\t' "$PL/sunflower_pruned.bim" > "$PL/sunflower_pruned.bim.tmp" \
  && mv "$PL/sunflower_pruned.bim.tmp" "$PL/sunflower_pruned.bim"

wc -l "$PL/sunflower.bim" "$PL/sunflower_pruned.bim" "$PL/sunflower.fam"
echo "DONE 02"
