#!/bin/bash
#SBATCH --job-name=sf04_adm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-120                       # KMAX(12) * NREP(10); see mapping below
#SBATCH --output=logs/sf04_adm_%A_%a.out
#SBATCH --error=logs/sf04_adm_%A_%a.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# task id -> (K, replicate)
i=$((SLURM_ARRAY_TASK_ID - 1))
K=$(( i / NREP + KMIN ))
REP=$(( i % NREP + 1 ))
SEED=$(( 1000 * K + REP ))

AD=$OUT/03_admixture/K${K}_rep${REP}
mkdir -p "$AD"; cd "$AD"

BED=$OUT/01_qc/sunflower_pruned
ln -sf "$BED".bed .; ln -sf "$BED".bim .; ln -sf "$BED".fam .

admixture --cv="$ADMIX_CV" -j"$THREADS" -s "$SEED" \
          sunflower_pruned.bed "$K" | tee "log_K${K}_rep${REP}.out"

echo "DONE 04 K=$K rep=$REP"
