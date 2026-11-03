#!/bin/bash
#SBATCH --job-name=sf17_depth_ld
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/sf17_depth_ld_%j.out
#SBATCH --error=logs/sf17_depth_ld_%j.err
set -euo pipefail
source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Read depth at the sites removed by the heterozygosity ceiling (Table S3f) and
# linkage disequilibrium corrected for kinship (Table S14d), both from the
# delivered matrices. The kinship matrix, the accession table, the ancestry
# proportions and the list of retained sites are read from $OUT/12_review/ld_inputs.
cd "$ANA"
python3 -c "import numpy, pandas" 2>/dev/null || python3 -m pip install numpy pandas
python3 -c "import sys, numpy, pandas; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'pandas', pandas.__version__)"
IN=$OUT/12_review/ld_inputs
python3 scripts/site_depth_by_filter.py "$RAW/Sunflower.MCR50.AlleleCounts" "$RAW/Sunflower.MCR50.snps.genotype.txt" "$OUT/01_qc/keep_samples.txt" "$OUT/09_tables"
python3 scripts/ld_kinship_corrected.py "$RAW/Sunflower.MCR50.snps.Imputed.genotypes.txt" "$IN/Data1_kinship.csv" "$IN/TableS1_accession_list.csv" "$IN/TableS6_admixture_Q.csv" "$IN/retained_sites.txt" "$OUT/09_tables"
echo "DONE 17"
