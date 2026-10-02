#!/usr/bin/env bash
# ============================================================================
# Sunflower tGBS - Genetic Diversity & Population Structure
# Global configuration. Every script sources this file.
#   usage:  source /path/to/analysis/config.sh
# ============================================================================

# ---- 1. Project roots (Swan: ydelen2, group dweikat) -----------------------
# /work        = Lustre scratch. Fast parallel I/O, PURGED after 6 months, not
#                backed up.  Run the jobs here.
# /mnt/nrdstor = NRDStor persistent research storage.  Keep the Freedom Markers
#                delivery and the final figures/tables here; do not run jobs
#                against it (metadata-heavy SLURM I/O is slow on NRDStor).
export PROJ=${PROJ:-/work/dweikat/ydelen2/sunflower_tgbs}
export ARCHIVE=${ARCHIVE:-/mnt/nrdstor/dweikat/ydelen2/sunflower_tgbs}
export RAW=$PROJ/raw                 # copy of tables/ from the Freedom Markers drive
export REFDIR=$PROJ/genome           # HanXRQr2.0-SUNRISE.fasta
export ANA=$PROJ/analysis            # this directory (scripts)
export OUT=$PROJ/results
export LOGS=$PROJ/logs
export TMPD=${SLURM_SUBMIT_DIR:-$PROJ}/tmp

mkdir -p "$OUT"/{00_vcf,01_qc,02_pca,03_admixture,04_diversity,05_ld,06_tree,07_core,08_figures} \
         "$LOGS" "$TMPD"

# ---- 2. Input files --------------------------------------------------------
export GENO_MCR50=$RAW/Sunflower.MCR50.snps.genotype.txt        # 276 samples x 247,008 SNP (0/1/2/3 codes)
export GENO_ALL=$RAW/Sunflower.all.snps.genotype.txt            # 287 samples x 503,188 SNP
export VCF_IMP=$RAW/Sunflower.MCR50.snps.Imputed.vcf            # Beagle 5.1 imputed, phased, GT only
export HMP_IMP=$RAW/Sunflower.MCR50.snps.Imputed.hmp.txt
export FASTA=$REFDIR/HanXRQr2.0-SUNRISE.fasta
export CHRMAP=$ANA/ref/chr_map.txt
export CHRLEN=$ANA/ref/chr_lengths.txt
export META=$ANA/metadata/sample_metadata.csv                   # built by R/01_build_metadata.R

# ---- 3. Filtering thresholds ---------------------------------------------
export MAX_SAMPLE_MISS=0.90    # matches the accession filter of the published papers
                               # on this panel (274 accessions retained)
export MAX_SITE_MISS=0.50      # no filter beyond the provider MCR50 call rate
export MIN_MAF=0.05            # minor allele frequency
export MAX_SITE_HET=0.55       # HWE ceiling: He cannot exceed 0.50 at a biallelic locus,
                               # so >55% het indicates a collapsed duplicate. A single
                               # marker set is used for every analysis.
export LD_PRUNE="50 10 0.2"    # plink --indep-pairwise window step r2  (for PCA/ADMIXTURE/tree)
export KING_DUP=0.354          # KING kinship above this = duplicate / same accession

# ---- 4. ADMIXTURE grid -----------------------------------------------------
export KMIN=1
export KMAX=12
export NREP=10                 # independent seeds per K
export ADMIX_CV=10             # cross-validation folds

# ---- 5. Compute ------------------------------------------------------------
export THREADS=${SLURM_CPUS_PER_TASK:-8}
export PLINK_MEM=${PLINK_MEM:-32000}   # MB given to plink --memory

# ---- 6. Software -----------------------------------------------------------
# Preferred: one conda env holding everything (see env/environment.yml).
# Fallback:  HCC modules. Run `module spider <name>` on Swan to get exact names.
if [[ -z "${SUNFLOWER_CONDA:-}" && -d /work/dweikat/ydelen2/envs/sunflower ]]; then
    export SUNFLOWER_CONDA=/work/dweikat/ydelen2/envs/sunflower
fi

load_tools () {
  # HCC's conda.sh reads PS1, which is fatal under `set -u`. Relax nounset for
  # the duration of activation and restore whatever the caller had, so sourcing
  # this file in an interactive shell does not leave nounset switched on.
  local _nounset_was_set=0
  case $- in *u*) _nounset_was_set=1 ;; esac
  set +u; export PS1="${PS1:-}"
  if [[ -n "${SUNFLOWER_CONDA:-}" ]]; then
      module load anaconda 2>/dev/null || module load mamba 2>/dev/null || true
      # shellcheck disable=SC1091
      source activate "$SUNFLOWER_CONDA"
  else
      module load bcftools   2>/dev/null || echo "WARN: bcftools module not found"
      module load vcftools   2>/dev/null || echo "WARN: vcftools module not found"
      module load plink      2>/dev/null || echo "WARN: plink module not found"
      module load admixture  2>/dev/null || echo "WARN: admixture module not found"
      module load samtools   2>/dev/null || true
      module load R          2>/dev/null || echo "WARN: R module not found"
  fi
  [[ $_nounset_was_set -eq 1 ]] && set -u
  return 0
}
