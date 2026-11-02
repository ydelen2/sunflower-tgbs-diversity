## Step 23 - inputs for the review checks ------------------------------------------
## Writes the list of accessions dropped for the ADMIXTURE run without close
## relatives (one member of every pair with a genomic relationship above 0.5,
## the member with the higher missing rate) and summarizes relatedness within
## and between the ancestry clusters from the same relationship matrix
## (Data file 1, plink2 --make-grm-list on the imputed filtered set).
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
options(width = 200)
RV <- file.path(OUT, "12_review"); dir.create(RV, showWarnings = FALSE, recursive = TRUE)
PL <- file.path(OUT, "01_qc")

find_tab <- function(pattern) {
  f <- list.files(TAB, pattern, full.names = TRUE)
  if (length(f) != 1L) stop("expected one file matching ", pattern, ", found ", length(f))
  f
}

gl  <- fread(file.path(PL, "grm.grm"), col.names = c("i", "j", "nsnp", "k"))
ids <- fread(file.path(PL, "grm.grm.id"), header = FALSE)$V2
qc  <- fread(find_tab("_accession_qc\\.csv$"))
asg <- fread(file.path(OUT, "03_admixture/group_assignment.txt"),
             header = FALSE, col.names = c("SAMPLE", "GROUP"))
miss <- qc$MISS_RATE[match(ids, qc$SAMPLE)]
grp  <- asg$GROUP[match(ids, asg$SAMPLE)]
stopifnot(!anyNA(miss), !anyNA(grp))

## ---- relatedness within and between clusters ---------------------------------
off <- gl[i != j]
off[, `:=`(G1 = grp[i], G2 = grp[j])]
within <- off[G1 == G2, .(N = length(unique(c(i, j))),
                          WITHIN_MEAN = mean(k), WITHIN_MEDIAN = median(k), WITHIN_MAX = max(k),
                          PAIRS = .N, PAIRS_ABOVE_0_25 = sum(k > 0.25),
                          PAIRS_ABOVE_0_5 = sum(k > 0.5), PAIRS_ABOVE_0_75 = sum(k > 0.75)),
              by = .(GROUP = G1)]
between <- off[G1 != G2, .(BETWEEN_MEAN = mean(k)), by = .(GROUP = G1)]
between <- rbind(between, off[G1 != G2, .(BETWEEN_MEAN = mean(k)), by = .(GROUP = G2)])[
  , .(BETWEEN_MEAN = mean(BETWEEN_MEAN)), by = GROUP]
rel <- merge(within, between, by = "GROUP")[order(GROUP)]
num <- setdiff(names(rel), c("GROUP", "N", "PAIRS", "PAIRS_ABOVE_0_25", "PAIRS_ABOVE_0_5", "PAIRS_ABOVE_0_75"))
rel[, (num) := lapply(.SD, round, 4), .SDcols = num]
print(rel)
fwrite(rel, file.path(TAB, "review_relatedness_by_cluster.csv"))

## ---- one accession dropped per closely related pair ---------------------------
pairs <- off[k > 0.5][order(-k)]
drop <- character(0)
for (r in seq_len(nrow(pairs))) {
  a <- ids[pairs$i[r]]; b <- ids[pairs$j[r]]
  if (a %in% drop || b %in% drop) next
  drop <- c(drop, if (miss[match(a, ids)] > miss[match(b, ids)]) a else b)
}
cat(sprintf("pairs with relationship > 0.5: %d | accessions dropped: %d\n", nrow(pairs), length(drop)))
print(table(grp[match(drop, ids)]))
writeLines(sort(drop), file.path(RV, "drop_related.txt"))
fwrite(data.table(SAMPLE = sort(drop), GROUP = grp[match(sort(drop), ids)]),
       file.path(TAB, "review_dropped_relatives.csv"))
cat("step 23 done\n")
