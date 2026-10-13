#!/bin/bash
#SBATCH --job-name=sf09_all
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/sf09_all_%j.out
#SBATCH --error=logs/sf09_all_%j.err
set -euo pipefail

source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Compares the permissive ALL-SNP delivery (503,188 sites, provider call rate
# >=20%) against the MCR50 set used for the main analyses (247,008 sites, call
# rate >=50%), on the same 274 accessions and with the same site filters apart
# from missingness. Answers whether the extra, lower-coverage sites change the
# diversity and LD estimates or only add noise.

AC=$OUT/10_allsnp
mkdir -p "$AC"
QC=$OUT/01_qc
cd "$AC"

# ---- 1. per-site statistics on every ALL-SNP chromosome site ----------------
# max-site-miss 0.80 is the provider's own floor, so nothing is removed here on
# missingness; the thresholds are applied afterwards from the stats file.
python3 "$ANA/scripts/genotype2vcf.py" \
    -i "$GENO_ALL" \
    --chr-map "$CHRMAP" --chrom-only \
    --keep "$QC/keep_samples.txt" \
    --max-site-miss 0.80 \
    --min-maf "$MIN_MAF" \
    --max-site-het "$MAX_SITE_HET" \
    --stats "$AC/site_stats_all.tsv" \
  | bgzip -@ "$THREADS" -c > allsnp.filtered.vcf.gz
tabix -f -p vcf allsnp.filtered.vcf.gz

python3 "$ANA/scripts/sample_stats.py" \
    -i "$GENO_ALL" --chr-map "$CHRMAP" --core-call 0.90 \
    --keep "$QC/keep_samples.txt" \
    -o "$AC/sample_stats_all.tsv"

# ---- 2. marker count and mean statistics at a range of thresholds -----------
for f in "$QC/site_stats_filtered.tsv" "$AC/site_stats_all.tsv"; do
  tag=$([[ $f == *site_stats_all* ]] && echo ALL || echo MCR50)
  awk -F'\t' -v T="$tag" -v A="$MIN_MAF" -v H="$MAX_SITE_HET" '
    NR>1 && $5>=A && $7<=H {
      for (i = 1; i <= 7; i++) {
        t = i * 0.1 + 0.1
        if ($4 <= t) { n[i]++; maf[i]+=$5; he[i]+=$6; ho[i]+=$7 }
      }
    }
    END {
      for (i = 1; i <= 7; i++) if (n[i] > 0)
        printf "%s\t%.1f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n", T, i*0.1+0.1, n[i],
               maf[i]/n[i], he[i]/n[i], ho[i]/n[i], 1-(ho[i]/n[i])/(he[i]/n[i])
    }' "$f"
done | { printf 'SET\tMAX_SITE_MISS\tN_SNP\tMAF\tHE\tHO\tFIS\n'; cat; } \
  > "$AC/threshold_comparison.tsv"
cat "$AC/threshold_comparison.tsv"
mkdir -p "$OUT/09_tables"
cp "$AC/threshold_comparison.tsv" "$OUT/09_tables/TableS16_callrate_threshold_comparison.tsv"

# ---- 3. LD on the ALL-SNP set ----------------------------------------------
if command -v PopLDdecay >/dev/null 2>&1; then
  PopLDdecay -InVCF allsnp.filtered.vcf.gz -OutStat allsnp.stat \
             -MaxDist 1000 -MAF "$MIN_MAF" -Miss 0.8
  {
    printf 'SET\tDISTANCE_CLASS\tR2\tPAIRS\n'
    for s in allsnp "$OUT/05_ld/all"; do
      [[ -e ${s}.stat.gz ]] || continue
      lab=$([[ $(basename "$s") == allsnp ]] && echo "Permissive set, missing rate up to 0.80" || echo "Analysis set")
      zcat "${s}.stat.gz" | awk -v L="$lab" 'NR>1 && $6>0 {
          k = ($1<=2000)?0 : ($1<=10000)?1 : ($1<=50000)?2 : ($1<=200000)?3 : 4
          s[k]+=$2*$6; n[k]+=$6 }
        END { split("0-2kb 2-10kb 10-50kb 50-200kb 200kb-1Mb", o, " ")
              for (k=0;k<=4;k++) if (n[k])
                printf "%s\t%s\t%.4f\t%d\n", L, o[k+1], s[k]/n[k], n[k] }'
    done
  } > "$OUT/09_tables/TableS16b_callrate_LD_comparison.tsv"
  cat "$OUT/09_tables/TableS16b_callrate_LD_comparison.tsv"
fi

echo "DONE 09"
