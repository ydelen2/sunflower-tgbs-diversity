#!/bin/bash
#SBATCH --job-name=sf03_pca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sf03_pca_%j.out
#SBATCH --error=logs/sf03_pca_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

PL=$OUT/01_qc
PC=$OUT/02_pca

# PCA on the LD-pruned imputed set.
plink2 --bfile "$PL/sunflower_pruned" --allow-extra-chr \
       --pca 20 allele-wts \
       --memory "$PLINK_MEM" --threads "$THREADS" \
       --out "$PC/pca"

# The same PCA on the non-imputed genotypes at the same markers; PLINK 2
# mean-imputes the missing calls. Used only to check the effect of imputation
# on the ordination (R/17 correlates the two sets of components).
plink2 --bfile "$PL/sunflower" --allow-extra-chr        --extract "$PL/prune.prune.in"        --pca 20 allele-wts        --memory "$PLINK_MEM" --threads "$THREADS"        --out "$PC/pca_nonimputed"

# GRM on the same pruned markers as the PCA. The eigenvalues sum to this
# matrix's trace, so it is the correct denominator for variance explained;
# the GRM in step 02 is built on the full marker set and does not match.
plink2 --bfile "$PL/sunflower_pruned" --allow-extra-chr        --make-grm-list --memory "$PLINK_MEM" --threads "$THREADS"        --out "$PC/pca_grm"

# variance explained
awk 'NR==FNR{s+=$1; next} {printf "PC%d\t%.4f\t%.4f\n", FNR, $1, 100*$1/s}' \
    "$PC/pca.eigenval" "$PC/pca.eigenval" > "$PC/pca.varexp.tsv"

echo "DONE 03"
