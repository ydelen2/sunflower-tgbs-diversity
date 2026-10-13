## Step 21 - repeat of the core subset searches ------------------------------------
## Core Hunter's search is stochastic and its seed cannot be fixed from R, so
## every search of step 16 is run a second time with the same settings. The
## archived selections in results/07_core are left untouched; the repeat is
## written to results/07_core_repeat and the two runs are compared.
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
suppressPackageStartupMessages({ library(corehunter); library(rJava) })

PL  <- file.path(OUT, "01_qc")
CR  <- file.path(OUT, "07_core")
CR2 <- file.path(OUT, "07_core_repeat"); dir.create(CR2, showWarnings = FALSE, recursive = TRUE)

d <- fread(file.path(PL, "sunflower_pruned_A.raw"))
ids <- d$IID
m <- as.matrix(d[, -(1:6)])
rownames(m) <- ids
m[is.na(m)] <- round(colMeans(m, na.rm = TRUE))[col(m)][is.na(m)]
gdata <- genotypes(data = m, format = "biparental")
asg <- fread(file.path(OUT, "03_admixture/group_assignment.txt"),
             header = FALSE, col.names = c("SAMPLE", "GROUP"))

objs <- list(objective("EN", "MR"), objective("SH"))
evaluate <- function(sel) {
  data.table(N = length(sel),
             EN_MR   = evaluateCore(sel, gdata, objective("EN", "MR")),
             SHANNON = evaluateCore(sel, gdata, objective("SH")),
             HE      = evaluateCore(sel, gdata, objective("HE")))
}
fractions <- c(0.10, 0.15, 0.20, 0.25, 0.30)

## whole-panel searches, as in step 16
whole <- rbindlist(lapply(fractions, function(fr) {
  set.seed(99)
  core <- sampleCore(gdata, size = fr, obj = objs, mode = "default")
  writeLines(core$sel, file.path(CR2, sprintf("core_%02d.txt", round(fr * 100))))
  cbind(data.table(FRACTION = fr, METHOD = "Core Hunter, whole panel"), evaluate(core$sel))
}))

## stratified searches, as in step 16
strat_core <- function(fr, min_per_group = 2L) {
  sel <- character(0)
  for (g in unique(asg$GROUP)) {
    ids_g <- intersect(asg[GROUP == g, SAMPLE], ids)
    if (length(ids_g) == 0) next
    k <- max(min_per_group, round(fr * length(ids_g)))
    k <- min(k, length(ids_g))
    if (k >= length(ids_g)) { sel <- c(sel, ids_g); next }
    gd <- genotypes(data = m[ids_g, , drop = FALSE], format = "biparental")
    set.seed(99)
    cg <- sampleCore(gd, size = k, obj = objs, mode = "default")
    sel <- c(sel, cg$sel)
  }
  sel
}
strat <- rbindlist(lapply(fractions, function(fr) {
  sel <- strat_core(fr)
  writeLines(sel, file.path(CR2, sprintf("core_stratified_%02d.txt", round(fr * 100))))
  cbind(data.table(FRACTION = fr, METHOD = "stratified"), evaluate(sel))
}))
run2 <- rbind(whole, strat)

## the archived run, re-evaluated from its selections so that both runs are
## scored by the same code
archived <- rbindlist(lapply(fractions, function(fr) rbind(
  cbind(data.table(FRACTION = fr, METHOD = "Core Hunter, whole panel"),
        evaluate(readLines(file.path(CR, sprintf("core_%02d.txt", round(fr * 100)))))),
  cbind(data.table(FRACTION = fr, METHOD = "stratified"),
        evaluate(readLines(file.path(CR, sprintf("core_stratified_%02d.txt", round(fr * 100)))))))))

shared <- function(fr, method) {
  stem <- if (method == "stratified") "core_stratified_%02d.txt" else "core_%02d.txt"
  a <- readLines(file.path(CR,  sprintf(stem, round(fr * 100))))
  b <- readLines(file.path(CR2, sprintf(stem, round(fr * 100))))
  length(intersect(a, b))
}
cmp <- merge(archived, run2, by = c("FRACTION", "METHOD"), suffixes = c("_RUN1", "_RUN2"))
cmp[, `:=`(EN_MR_DIFF = EN_MR_RUN2 - EN_MR_RUN1,
           SHANNON_DIFF = SHANNON_RUN2 - SHANNON_RUN1,
           HE_DIFF = HE_RUN2 - HE_RUN1,
           SHARED_ACCESSIONS = mapply(shared, FRACTION, METHOD))]
cmp[, SHARED_PCT := round(100 * SHARED_ACCESSIONS / N_RUN1, 1)]
cmp <- cmp[order(METHOD, FRACTION),
           .(FRACTION, METHOD, N = N_RUN1, EN_MR_RUN1, EN_MR_RUN2, EN_MR_DIFF,
             SHANNON_RUN1, SHANNON_RUN2, SHANNON_DIFF, HE_RUN1, HE_RUN2, HE_DIFF,
             SHARED_ACCESSIONS, SHARED_PCT)]
fwrite(cmp, file.path(TAB, "TableS15c_core_repeat.csv"))
options(width = 200)
print(cmp[, .(FRACTION, METHOD, N, SHANNON_RUN1, SHANNON_RUN2, SHANNON_DIFF, HE_DIFF, EN_MR_DIFF, SHARED_ACCESSIONS, SHARED_PCT)], digits = 5)
cat(sprintf("largest absolute difference: Shannon %.5f, He %.5f, EN-MR %.5f\n",
            max(abs(cmp$SHANNON_DIFF)), max(abs(cmp$HE_DIFF)), max(abs(cmp$EN_MR_DIFF))))
cat("step 21 done\n")
