## Step 17 - the statistical tests quoted in the manuscript text ----------------
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))

QC  <- file.path(OUT, "01_qc")
pcs <- fread(file.path(TAB, "Data4_pca_scores.csv"), colClasses = "character")
roh <- fread(file.path(TAB, "Data7_ROH.csv"), colClasses = "character")
qc  <- fread(file.path(QC, "sample_stats_kept.tsv"))

num <- function(x) suppressWarnings(as.numeric(x))
pcs[, `:=`(PC1 = num(PC1), OIL_PCT = num(OIL_PCT))]
roh[, F_ROH := num(F_ROH)]

res <- list()
add <- function(test, term, statistic, df1, df2, p, n, note = "") {
  res[[length(res) + 1L]] <<- data.table(
    TEST = test, TERM = term, STATISTIC = statistic,
    DF1 = df1, DF2 = df2, P_VALUE = p, N = n, NOTE = note)
}

## 1. seed oil content by ancestry cluster
oil_grp <- merge(pcs[, .(SAMPLE, OIL_PCT)],
                 roh[, .(SAMPLE, GROUP)], by = "SAMPLE")[!is.na(OIL_PCT)]
a1 <- anova(lm(OIL_PCT ~ factor(GROUP), data = oil_grp))
add("Seed oil content by cluster", "Cluster", round(a1[1, "F value"], 2),
    a1[1, "Df"], a1[2, "Df"], signif(a1[1, "Pr(>F)"], 3), nrow(oil_grp),
    "one-way analysis of variance")

## 2. PC1 by curator-assigned market type
ty <- pcs[TYPE %chin% c("oilseed", "confection")]
if (ty[, uniqueN(TYPE)] == 2L) {
  t1 <- t.test(PC1 ~ TYPE, data = ty)
  add("PC1 by market type", "Oilseed versus confection", round(t1$statistic, 2),
      round(t1$parameter, 1), NA_real_, signif(t1$p.value, 3), nrow(ty),
      "Welch two-sample t test")
}

## 3. PC1 against seed oil content
c1 <- pcs[!is.na(OIL_PCT), cor.test(PC1, OIL_PCT)]
add("PC1 against seed oil content", "Pearson correlation", round(c1$estimate, 3),
    c1$parameter, NA_real_, signif(c1$p.value, 3), pcs[!is.na(OIL_PCT), .N],
    "two-sided test of zero correlation")

## 4. PC1 by heterotic group
hg <- pcs[HETEROTIC_GROUP %chin% c("R", "B")]
if (hg[, uniqueN(HETEROTIC_GROUP)] == 2L) {
  t2 <- t.test(PC1 ~ HETEROTIC_GROUP, data = hg)
  add("PC1 by heterotic group", "Maintainer versus restorer", round(t2$statistic, 2),
      round(t2$parameter, 1), NA_real_, signif(t2$p.value, 3), nrow(hg),
      "Welch two-sample t test")
}

## 5. genomic inbreeding against per-accession missing rate, and the cluster
## effect adjusted for it
rq <- merge(roh[, .(SAMPLE, GROUP, F_ROH)], qc[, .(SAMPLE, MISS_RATE)], by = "SAMPLE")
c2 <- rq[, cor.test(F_ROH, MISS_RATE)]
add("Genomic inbreeding against missing rate", "Pearson correlation",
    round(c2$estimate, 3), c2$parameter, NA_real_, signif(c2$p.value, 3), nrow(rq),
    "two-sided test of zero correlation")

m1 <- lm(F_ROH ~ MISS_RATE, data = rq)
s1 <- summary(m1)$coefficients
add("Genomic inbreeding against missing rate", "Regression slope",
    round(s1["MISS_RATE", "Estimate"], 3), 1, m1$df.residual,
    signif(s1["MISS_RATE", "Pr(>|t|)"], 3), nrow(rq), "simple linear regression")

m2 <- lm(F_ROH ~ MISS_RATE + factor(GROUP), data = rq)
a2 <- anova(m2)
add("Genomic inbreeding, cluster adjusted for missing rate", "Cluster",
    round(a2["factor(GROUP)", "F value"], 2), a2["factor(GROUP)", "Df"],
    m2$df.residual, signif(a2["factor(GROUP)", "Pr(>F)"], 3), nrow(rq),
    "analysis of covariance, missing rate as covariate")

s2 <- summary(m2)$coefficients
add("Genomic inbreeding, cluster adjusted for missing rate", "Missing rate slope",
    round(s2["MISS_RATE", "Estimate"], 3), 1, m2$df.residual,
    signif(s2["MISS_RATE", "Pr(>|t|)"], 3), nrow(rq),
    "slope within clusters, from the same model")

## cluster means adjusted to the panel mean missing rate
adj <- data.table(GROUP = sort(unique(rq$GROUP)))
adj[, MISS_RATE := mean(rq$MISS_RATE)]
adj[, F_ROH_ADJUSTED := round(predict(m2, newdata = adj), 3)]
adj[, F_ROH_OBSERVED := rq[, round(mean(F_ROH), 3), by = GROUP][match(adj$GROUP, GROUP), V1]]
setorder(adj, -F_ROH_ADJUSTED)
print(adj[, .(GROUP, F_ROH_OBSERVED, F_ROH_ADJUSTED)])

out <- rbindlist(res)
print(out)
fwrite(out, file.path(TAB, "TableS10_manuscript_statistics.csv"))
fwrite(adj[, .(GROUP, F_ROH_OBSERVED, F_ROH_ADJUSTED)],
       file.path(TAB, "TableS17_FROH_adjusted_means.csv"))

## sensitivity of the assignment to the ancestry threshold
qm <- fread(file.path(TAB, "TableS5_admixture_Q.csv"), colClasses = "character")
qm[, MAXQ := num(MAXQ)]
qm <- merge(qm, roh[, .(SAMPLE, HETEROTIC_GROUP)], by = "SAMPLE", all.x = TRUE)
thr <- rbindlist(lapply(c(0.50, 0.60, 0.70, 0.80, 0.90), function(t) {
  g <- fifelse(qm$MAXQ >= t, qm$CLUSTER, "admixed")
  cnt <- table(factor(g, levels = c("Q1", "Q2", "Q3", "Q4", "Q5", "admixed")))
  data.table(THRESHOLD = t,
             N_Q1 = cnt[["Q1"]], N_Q2 = cnt[["Q2"]], N_Q3 = cnt[["Q3"]],
             N_Q4 = cnt[["Q4"]], N_Q5 = cnt[["Q5"]], N_ADMIXED = cnt[["admixed"]],
             Q2_MAINTAINER = sum(g == "Q2" & qm$HETEROTIC_GROUP == "B", na.rm = TRUE),
             Q2_RESTORER   = sum(g == "Q2" & qm$HETEROTIC_GROUP == "R", na.rm = TRUE),
             Q5_MAINTAINER = sum(g == "Q5" & qm$HETEROTIC_GROUP == "B", na.rm = TRUE),
             Q5_RESTORER   = sum(g == "Q5" & qm$HETEROTIC_GROUP == "R", na.rm = TRUE))
}))
print(thr)
fwrite(thr, file.path(TAB, "TableS6_assignment_threshold.csv"))

## effect of imputation on the ordination
PC <- file.path(OUT, "02_pca")
f_imp <- file.path(PC, "pca.eigenvec"); f_non <- file.path(PC, "pca_nonimputed.eigenvec")
if (file.exists(f_imp) && file.exists(f_non)) {
  a <- fread(f_imp); b <- fread(f_non)
  setnames(a, 1:2, c("FID", "SAMPLE")); setnames(b, 1:2, c("FID", "SAMPLE"))
  ab <- merge(a, b, by = "SAMPLE", suffixes = c("_imp", "_non"))
  pcc <- rbindlist(lapply(1:5, function(k) data.table(
    COMPONENT = paste0("PC", k),
    CORRELATION = round(abs(cor(ab[[paste0("PC", k, "_imp")]], ab[[paste0("PC", k, "_non")]])), 4),
    N = nrow(ab))))
  print(pcc)
  fwrite(pcc, file.path(TAB, "pca_imputation_check.csv"))
}

## accession list: passport data, descriptors and assignments in one table
qm2 <- fread(file.path(TAB, "TableS5_admixture_Q.csv"), colClasses = "character")
dap <- fread(file.path(TAB, "Data3_dapc_coordinates.csv"), colClasses = "character")
k3  <- fread(file.path(TAB, "Data6_K3_assignment.csv"), colClasses = "character")
acc <- pcs[, .(SAMPLE, ACCESSION, NAME, COUNTRY, REGION, STATUS, ORIGIN_EVENT,
               MARKET_TYPE = TYPE, HETEROTIC_GROUP, OIL_PCT, YEAR)]
acc[, HETEROTIC_GROUP := fcase(HETEROTIC_GROUP == "B", "maintainer",
                               HETEROTIC_GROUP == "R", "restorer",
                               default = "")]
acc <- merge(acc, qm2[, .(SAMPLE, CLUSTER = GROUP, MAX_ANCESTRY = MAXQ)], by = "SAMPLE", all.x = TRUE)
acc <- merge(acc, dap[, .(SAMPLE, DAPC_GROUP = DAPC)], by = "SAMPLE", all.x = TRUE)
acc <- merge(acc, k3[, .(SAMPLE, K3_GROUP)], by = "SAMPLE", all.x = TRUE)
acc <- merge(acc, qc[, .(SAMPLE, CALL_RATE = round(CALL_RATE, 4), MISS_RATE = round(MISS_RATE, 4),
                         HO = round(HO, 4))],
             by = "SAMPLE", all.x = TRUE)
acc <- merge(acc, roh[, .(SAMPLE, F_ROH = round(F_ROH, 4))], by = "SAMPLE", all.x = TRUE)
for (cl in names(acc)) if (is.character(acc[[cl]])) acc[is.na(get(cl)), (cl) := ""]
setorder(acc, CLUSTER, SAMPLE)
stopifnot(nrow(acc) == 274L)
fwrite(acc, file.path(TAB, "TableS1_accession_list.csv"))

cat("step 17 done\n")
