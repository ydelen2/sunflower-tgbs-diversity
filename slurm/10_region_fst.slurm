#!/bin/bash
#SBATCH --job-name=sf10_reg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sf10_reg_%j.out
#SBATCH --error=logs/sf10_reg_%j.err
set -euo pipefail

source "${ANA_CFG:-${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
load_tools

# Differentiation among geographic regions, measured directly rather than as a
# variance component conditioned on the ancestry clusters. Regions are tested
# across the whole panel and, separately, within the accessions whose GRIN
# origin field records a collection event, where the country is the place the
# material came from rather than the place it was bred.

RG=$OUT/11_region
mkdir -p "$RG/pops"
cd "$RG"

VCF=$OUT/00_vcf/sunflower.filtered.vcf.gz
KEEP=$OUT/01_qc/keep_samples.txt
META=$ANA/metadata/sample_metadata.csv
MINN=10

# region and origin-event columns, restricted to the analysis accessions
python3 - "$META" "$KEEP" > accession_region.tsv <<'PY'
import csv, sys
keep = {l.strip() for l in open(sys.argv[2]) if l.strip()}
for r in csv.DictReader(open(sys.argv[1], encoding="utf-8")):
    if r["SAMPLE"] in keep:
        reg = (r.get("REGION") or "").strip().replace(" ", "_").replace("/", "-")
        ev = (r.get("ORIGIN_EVENT") or "").strip()
        if reg and reg != "Unknown":
            print("\t".join([r["SAMPLE"], reg, ev]))
PY
wc -l < accession_region.tsv

run_set () {                       # $1 = label, $2 = awk filter on the event column
  local tag=$1 filt=$2
  rm -f pops/*."${tag}".ind
  awk -F'\t' -v F="$filt" -v T="$tag" \
      '$3 ~ F { print $1 > ("pops/" $2 "." T ".ind") }' \
      accession_region.tsv
  : > "pairwise_fst_${tag}.tsv"
  printf 'REGION1\tREGION2\tN1\tN2\tMEAN_FST\tWEIGHTED_FST\n' > "pairwise_fst_${tag}.tsv"
  mapfile -t R < <(for f in pops/*."${tag}".ind; do
                     [[ -e $f ]] || continue
                     [[ $(wc -l < "$f") -ge $MINN ]] && basename "$f" ".${tag}.ind"
                   done)
  echo "== $tag: ${#R[@]} regions with at least $MINN accessions: ${R[*]}"
  for ((i=0; i<${#R[@]}; i++)); do
    for ((j=i+1; j<${#R[@]}; j++)); do
      a=${R[$i]}; b=${R[$j]}
      vcftools --gzvcf "$VCF" \
               --weir-fst-pop "pops/${a}.${tag}.ind" \
               --weir-fst-pop "pops/${b}.${tag}.ind" \
               --out "fst_${tag}_${a}_vs_${b}" 2> "fst_${tag}_${a}_vs_${b}.log"
      m=$(awk '/Weir and Cockerham mean Fst estimate/{v=$NF} END{print v}'     "fst_${tag}_${a}_vs_${b}.log")
      w=$(awk '/Weir and Cockerham weighted Fst estimate/{v=$NF} END{print v}' "fst_${tag}_${a}_vs_${b}.log")
      printf '%s\t%s\t%d\t%d\t%s\t%s\n' "$a" "$b" \
             "$(wc -l < "pops/${a}.${tag}.ind")" "$(wc -l < "pops/${b}.${tag}.ind")" \
             "${m:-NA}" "${w:-NA}" >> "pairwise_fst_${tag}.tsv"
    done
  done
  column -t "pairwise_fst_${tag}.tsv"
}

run_set all       '.'
run_set collected 'Collected'

mkdir -p "$OUT/09_tables"
cp pairwise_fst_all.tsv       "$OUT/09_tables/TableS11_region_fst.tsv"
cp pairwise_fst_collected.tsv "$OUT/09_tables/TableS11b_region_fst_collected.tsv"

echo "DONE 10"
