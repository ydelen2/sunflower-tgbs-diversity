#!/bin/bash
#SBATCH --job-name=sf11_ldsens
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sf11_ldsens_%j.out
#SBATCH --error=logs/sf11_ldsens_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Two sensitivity runs for the whole-panel LD curve, summarized by R/15:
#   sens_nonimputed  the non-imputed filtered set, same sites and accessions
#   sens_het020      the imputed set restricted to sites with a heterozygous
#                    call rate of at most 0.20, a stricter paralog filter
LD=$OUT/05_ld
VC=$OUT/00_vcf
QC=$OUT/01_qc
MAXDIST=1000
cd "$LD"

command -v PopLDdecay >/dev/null 2>&1 || { echo "PopLDdecay not found"; exit 1; }

PopLDdecay -InVCF "$VC/sunflower.filtered.vcf.gz" -OutStat sens_nonimputed.stat \
           -MaxDist $MAXDIST -MAF "$MIN_MAF" -Miss "$MAX_SITE_MISS"

# site_stats_filtered.tsv: CHROM POS N_CALLED MISS MAF HE HO FIS PASS REASON
awk -F'\t' 'NR>1 && $9==1 && $7<=0.20 {print $1"\t"$2}' "$QC/site_stats_filtered.tsv" \
  > sens_het020.sites.tsv
printf 'sites with heterozygous call rate <= 0.20: %s\n' "$(wc -l < sens_het020.sites.tsv)"
bcftools view -T sens_het020.sites.tsv "$VC/sunflower.filtered.imputed.vcf.gz" \
  -Oz -o sens_het020.vcf.gz --threads "$THREADS"
tabix -f -p vcf sens_het020.vcf.gz
PopLDdecay -InVCF sens_het020.vcf.gz -OutStat sens_het020.stat \
           -MaxDist $MAXDIST -MAF "$MIN_MAF" -Miss 0.2

echo "DONE 11"
