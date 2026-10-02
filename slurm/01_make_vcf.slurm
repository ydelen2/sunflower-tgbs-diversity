#!/bin/bash
#SBATCH --job-name=sf01_vcf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/sf01_vcf_%j.out
#SBATCH --error=logs/sf01_vcf_%j.err
##SBATCH --partition=batch        # uncomment/edit after checking `sinfo -s` on Swan
set -euo pipefail

source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

cd "$OUT/00_vcf"

# ---------------------------------------------------------------------------
# 1. Rebuild a clean VCF from the genotype matrix (chromosomes 1-17 only).
#    No site filters yet; accession QC is judged on the full matrix.
# ---------------------------------------------------------------------------
python3 "$ANA/scripts/genotype2vcf.py" \
    -i "$GENO_MCR50" \
    --chr-map "$CHRMAP" --chrom-only \
    --min-called 2 \
    --stats "$OUT/01_qc/site_stats_raw.tsv" \
  | bgzip -@ "$THREADS" -c > mcr50.chrom.vcf.gz
tabix -f -p vcf mcr50.chrom.vcf.gz

# ---------------------------------------------------------------------------
# 2. Per-accession QC table (raw Ho + Ho on a common well-covered site subset)
# ---------------------------------------------------------------------------
python3 "$ANA/scripts/sample_stats.py" \
    -i "$GENO_MCR50" --chr-map "$CHRMAP" --core-call 0.90 \
    -o "$OUT/01_qc/sample_stats.tsv"

# ---------------------------------------------------------------------------
# 3. The Beagle-imputed set: rename contigs, keep chromosomes only.
#    Done with awk because the delivered VCF carries no ##contig header lines,
#    which bcftools annotate refuses to rename.  awk also injects the missing
#    contig headers so that downstream bcftools calls accept the file.
# ---------------------------------------------------------------------------
awk -v OFS='\t' '
  NR==FNR       { map[$1]=$2; order[++k]=$2; next }
  /^##/         { print; next }
  /^#CHROM/     { for (i=1;i<=k;i++) printf "##contig=<ID=%s>\n", order[i]; print; next }
  ($1 in map)   { $1=map[$1]; $3=$1"_"$2; print }' \
  "$CHRMAP" "$VCF_IMP" \
  | bgzip -@ "$THREADS" -c > mcr50.imputed.chrom.vcf.gz
tabix -f -p vcf mcr50.imputed.chrom.vcf.gz

# ---------------------------------------------------------------------------
# 4. Per-accession QC on the permissive ALL-SNP set (287 libraries), reported
#    for comparison with the MCR50 set
# ---------------------------------------------------------------------------
if [[ -s "$GENO_ALL" ]]; then
  python3 "$ANA/scripts/sample_stats.py" \
      -i "$GENO_ALL" --chr-map "$CHRMAP" --core-call 0.90 \
      -o "$OUT/01_qc/sample_stats_allsnps_287.tsv"
fi

echo "DONE 01"
