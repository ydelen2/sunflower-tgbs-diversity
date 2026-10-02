#!/usr/bin/env bash
# Fetch the HanXRQr2.0-SUNRISE (GCF_002127325.2) annotation from NCBI and reduce
# it to a few derived files.  Run once, on a login node
# (Swan login nodes have outbound internet; compute nodes may not).
#
#   bash scripts/get_annotation.sh
set -euo pipefail
source "${ANA_CFG:-$PWD/config.sh}"

AN=$REFDIR/annotation
mkdir -p "$AN"; cd "$AN"

BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/127/325/GCF_002127325.2_HanXRQr2.0-SUNRISE
GFF=GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.gff.gz

if [[ ! -s $GFF ]]; then
  echo "downloading annotation from NCBI ..."
  curl -fL -O "$BASE/$GFF"
  curl -fL -O "$BASE/md5checksums.txt" || true
fi
[[ -s md5checksums.txt ]] && grep -F "$GFF" md5checksums.txt | md5sum -c - || \
  echo "NOTE: checksum file unavailable, skipping verification"

# ---------------------------------------------------------------------------
# genes.bed  : chrom(1-17)  start0  end  gene_id  .  strand  description
# ---------------------------------------------------------------------------
zcat "$GFF" \
  | awk -v OFS='\t' -F'\t' '$3=="gene"{
        id="."; nm=".";
        n=split($9,a,";");
        for(i=1;i<=n;i++){
          if(a[i]~/^ID=/)   {id=substr(a[i],4)}
          if(a[i]~/^gene=/) {nm=substr(a[i],6)}
        }
        print $1, $4-1, $5, nm, id, $7 }' \
  | awk -v OFS='\t' 'NR==FNR{m[$1]=$2; next} ($1 in m){$1=m[$1]; print}' "$CHRMAP" - \
  | sort -k1,1n -k2,2n > genes.bed

# ---------------------------------------------------------------------------
# gene_density_1Mb.tsv : chrom  win_start  n_genes  coding_bp
#   Used to locate pericentromeric regions without a genetic map: in sunflower,
#   as in most large plant genomes, gene density collapses across the
#   pericentromere.
# ---------------------------------------------------------------------------
awk -v OFS='\t' -v W=1000000 '{
      w=int($2/W); n[$1"\t"w]++; bp[$1"\t"w]+=$3-$2 }
    END{ print "CHROM","WIN_START","N_GENES","GENIC_BP";
         for(k in n){ split(k,a,"\t"); print a[1], a[2]*W, n[k], bp[k] } }' genes.bed \
  | (read -r h; echo "$h"; sort -k1,1n -k2,2n) > gene_density_1Mb.tsv

wc -l genes.bed gene_density_1Mb.tsv
echo "annotation ready in $AN"
