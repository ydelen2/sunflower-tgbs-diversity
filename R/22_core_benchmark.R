## Step 22 - benchmark of the core subsets against simple sampling schemes -------
## Random subsets, proportional stratified random subsets (same allocation as
## the stratified Core Hunter scheme) and the two archived Core Hunter subsets,
## all scored with the same Core Hunter objectives on the pruned marker set,
## plus the share of the panel's alleles retained. 100 replicates per scheme.
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
suppressPackageStartupMessages({ library(corehunter); library(rJava) })

PL <- file.path(OUT, "01_qc")
CR <- file.path(OUT, "07_core")
d <- fread(file.path(PL, "sunflower_pruned_A.raw"))
ids <- d$IID
m <- as.matrix(d[, -(1:6)])
rownames(m) <- ids
m[is.na(m)] <- round(colMeans(m, na.rm = TRUE))[col(m)][is.na(m)]
gdata <- genotypes(data = m, format = "biparental")
asg <- fread(file.path(OUT, "03_admixture/group_assignment.txt"),
             header = FALSE, col.names = c("SAMPLE", "GROUP"))
asg <- asg[SAMPLE %in% ids]

allele_cov <- function(sel) {
  sub <- m[sel, , drop = FALSE]
  present <- (colSums(sub > 0) > 0) + (colSums(sub < 2) > 0)
  sum(present) / (2 * ncol(m))
}
score <- function(sel, scheme, rep) {
  data.table(SCHEME = scheme, REP = rep, N = length(sel),
             EN_MR = evaluateCore(sel, gdata, objective("EN", "MR")),
             SHANNON = evaluateCore(sel, gdata, objective("SH")),
             HE = evaluateCore(sel, gdata, objective("HE")),
             CV = evaluateCore(sel, gdata, objective("CV")),      # Core Hunter allele coverage
             ALLELE_COVERAGE = allele_cov(sel))
}

fr <- 0.20
n_total <- round(fr * length(ids))
alloc <- asg[, .(k = pmin(.N, pmax(2L, round(fr * .N)))), by = GROUP]
set.seed(2026)
res <- list()
for (rep in 1:100) {
  res[[length(res) + 1]] <- score(sample(ids, n_total), "random", rep)
  sel <- unlist(lapply(seq_len(nrow(alloc)), function(i)
    sample(asg[GROUP == alloc$GROUP[i], SAMPLE], alloc$k[i])))
  res[[length(res) + 1]] <- score(sel, "proportional stratified random", rep)
}
res[[length(res) + 1]] <- score(readLines(file.path(CR, "core_20.txt")), "Core Hunter, whole panel", 0)
res[[length(res) + 1]] <- score(readLines(file.path(CR, "core_stratified_20.txt")), "Core Hunter, stratified", 0)
res[[length(res) + 1]] <- score(ids, "full collection", 0)
out <- rbindlist(res)
fwrite(out, file.path(TAB, "review_core_benchmark_replicates.csv"))
summ <- out[, .(N = N[1], REPLICATES = .N,
                EN_MR = mean(EN_MR), EN_MR_SD = sd(EN_MR),
                SHANNON = mean(SHANNON), SHANNON_SD = sd(SHANNON),
                HE = mean(HE), HE_SD = sd(HE), CV = mean(CV), CV_SD = sd(CV),
                ALLELE_COVERAGE = mean(ALLELE_COVERAGE), ALLELE_COVERAGE_SD = sd(ALLELE_COVERAGE)),
            by = SCHEME]
fwrite(summ, file.path(TAB, "review_core_benchmark.csv"))
options(width = 200)
print(summ, digits = 5)
## where the two Core Hunter subsets fall in the random distributions
for (s in c("random", "proportional stratified random")) {
  r <- out[SCHEME == s]
  for (cs in c("Core Hunter, whole panel", "Core Hunter, stratified")) {
    v <- out[SCHEME == cs]
    cat(sprintf("%-26s against %-32s: Shannon higher than %3d%% of replicates, He %3d%%, EN-MR %3d%%, allele coverage %3d%%\n",
                cs, s, round(100 * mean(r$SHANNON < v$SHANNON)), round(100 * mean(r$HE < v$HE)),
                round(100 * mean(r$EN_MR < v$EN_MR)), round(100 * mean(r$ALLELE_COVERAGE < v$ALLELE_COVERAGE))))
  }
}
cat("step 22 done\n")
