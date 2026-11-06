## Step 24 - analyses added for the first round of review ------------------------
## Each section reads saved results (and the files written by
## slurm/13_review_checks.slurm into results/12_review) and writes a
## review_*.csv table into results/09_tables. Sections are independent; a
## missing input skips its section with a message.
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
options(width = 200, datatable.print.nrows = 300)
RV <- file.path(OUT, "12_review")
DV <- file.path(OUT, "04_diversity")
AD <- file.path(OUT, "03_admixture")
PL <- file.path(OUT, "01_qc")
LD <- file.path(OUT, "05_ld")
CL <- paste0("Q", 1:5)

find_tab <- function(pattern) {
  f <- list.files(TAB, pattern, full.names = TRUE)
  if (length(f) != 1L) stop("expected one file matching ", pattern, ", found ", length(f))
  f
}
section <- function(s) cat("\n==", s, "\n")
rnd <- function(dt, digits = 4) {
  num <- names(dt)[vapply(dt, is.double, logical(1))]
  dt[, (num) := lapply(.SD, round, digits), .SDcols = num]
  dt[]
}

asg <- fread(file.path(AD, "group_assignment.txt"), header = FALSE, col.names = c("SAMPLE", "GROUP"))
acc <- fread(find_tab("_accession_list\\.csv$"), colClasses = "character", na.strings = c("", "NA"))
for (cl in intersect(c("OIL_PCT", "MISS_RATE", "F_ROH", "MAX_ANCESTRY", "CALL_RATE", "HO"), names(acc)))
  acc[, (cl) := as.numeric(get(cl))]
acc[, GROUP := asg$GROUP[match(SAMPLE, asg$SAMPLE)]]
stopifnot(!anyNA(acc$GROUP))

## ===========================================================================
## A. Genotype matrix of the non-imputed filtered set (VCFtools --012)
## ===========================================================================
section("A. genotype matrix")
f012 <- file.path(RV, "all.012")
if (!file.exists(f012)) stop("run slurm/13_review_checks.slurm first: ", f012, " is missing")
indv <- readLines(paste0(f012, ".indv"))
pos  <- fread(paste0(f012, ".pos"), header = FALSE, col.names = c("CHROM", "POS"))
L    <- nrow(pos)
G <- matrix(scan(f012, what = integer(), sep = "\t", quiet = TRUE), nrow = length(indv), byrow = TRUE)
G <- G[, -1, drop = FALSE]                 # the first field of each row is the row index
stopifnot(ncol(G) == L)
G[G == -1L] <- NA_integer_
rownames(G) <- indv
stopifnot(all(acc$SAMPLE %in% indv))
cat(sprintf("%d accessions x %d sites; %.1f%% missing\n", nrow(G), L, 100 * mean(is.na(G))))

site_stats <- function(ids) {
  sub <- G[ids, , drop = FALSE]
  n   <- colSums(!is.na(sub))
  het <- colSums(sub == 1L, na.rm = TRUE)
  alt <- colSums(sub, na.rm = TRUE)
  p   <- ifelse(n > 0, alt / (2 * n), NA_real_)
  ho  <- ifelse(n > 0, het / n, NA_real_)
  list(n = n, p = p, ho = ho)
}

## ===========================================================================
## B. Diversity indices on one site set per group
##    He is the unbiased within-population gene diversity of Nei and Chesser
##    (1983), Ho the mean over sites of the heterozygous share of called
##    genotypes, F_IS = 1 - Ho/He, all over the same sites. The published Table 1
##    took Ho from VCFtools --het, which skips sites fixed within the group, and
##    He from VCFtools --freq over all called sites; the two columns were on
##    different site sets. The per-accession VCFtools-style value is kept here
##    to show that it reproduces the published column.
## ===========================================================================
section("B. diversity indices on a unified site set")
div_row <- function(ids, label) {
  s <- site_stats(ids); n <- s$n; p <- s$p; ho <- s$ho
  he_b <- 2 * p * (1 - p)
  hs   <- ifelse(n > 1, n / (n - 1) * (he_b - ho / (2 * n)), NA_real_)
  maf  <- pmin(p, 1 - p)
  pic  <- 1 - (p^2 + (1 - p)^2) - 2 * p^2 * (1 - p)^2
  called <- n > 0
  polys  <- called & p > 0 & p < 1
  sub_all  <- G[ids, called, drop = FALSE]
  sub_poly <- G[ids, polys, drop = FALSE]
  ho_ind_all  <- mean(rowSums(sub_all == 1L, na.rm = TRUE) / rowSums(!is.na(sub_all)))
  ho_ind_poly <- mean(rowSums(sub_poly == 1L, na.rm = TRUE) / rowSums(!is.na(sub_poly)))
  rbindlist(lapply(list(c("all called sites", "called"), c("sites polymorphic within the group", "polys")), function(w) {
    sel <- if (w[2] == "called") called else polys
    data.table(GROUP = label, N = length(ids), SITES = w[1], N_SITES = sum(sel),
               PCT_POLYMORPHIC = 100 * sum(polys) / sum(called),
               MAF = mean(maf[sel]), HE_2PQ = mean(he_b[sel]), HE = mean(hs[sel], na.rm = TRUE),
               HO = mean(ho[sel]), FIS = 1 - mean(ho[sel]) / mean(hs[sel], na.rm = TRUE),
               PIC = mean(pic[sel]),
               HO_PER_ACCESSION_SAME_SITES = if (w[2] == "called") ho_ind_all else ho_ind_poly)
  }))
}
groups <- c(as.list(setNames(CL, CL)), list(admixed = "admixed"))
uni <- rbindlist(c(
  lapply(names(groups), function(g) div_row(acc[GROUP == g, SAMPLE], g)),
  list(div_row(acc$SAMPLE, "All accessions"),
       div_row(acc[ORIGIN_EVENT == "Collected", SAMPLE], "Collected material"),
       div_row(acc[ORIGIN_EVENT == "Developed", SAMPLE], "Developed material"))))
## the whole panel under the stricter paralog ceiling of 0.20
s_all <- site_stats(acc$SAMPLE)
strict <- s_all$ho <= 0.20
hs_all <- s_all$n / (s_all$n - 1) * (2 * s_all$p * (1 - s_all$p) - s_all$ho / (2 * s_all$n))
for (w in list(c("0.55 (analysis set)", "all"), c("0.20", "strict"))) {
  sel <- if (w[2] == "all") rep(TRUE, L) else strict
  cat(sprintf("heterozygous call rate ceiling %s: %d sites, Ho %.4f, He %.4f, F_IS %.4f\n", w[1], sum(sel),
              mean(s_all$ho[sel]), mean(hs_all[sel]), 1 - mean(s_all$ho[sel]) / mean(hs_all[sel])))
}
uni <- rbind(uni, data.table(GROUP = "All accessions, ceiling 0.20", N = nrow(acc), SITES = "sites with heterozygous call rate <= 0.20",
                             N_SITES = sum(strict), PCT_POLYMORPHIC = 100, MAF = mean(pmin(s_all$p, 1 - s_all$p)[strict]),
                             HE_2PQ = mean((2 * s_all$p * (1 - s_all$p))[strict]), HE = mean(hs_all[strict]), HO = mean(s_all$ho[strict]),
                             FIS = 1 - mean(s_all$ho[strict]) / mean(hs_all[strict]), PIC = NA_real_, HO_PER_ACCESSION_SAME_SITES = NA_real_))
print(rnd(copy(uni)))
fwrite(rnd(copy(uni)), file.path(TAB, "review_diversity_unified.csv"))

t1 <- fread(find_tab("_cluster_diversity_and_identity\\.csv$"))
cat("\npublished Table 1 Ho against the per-accession VCFtools-style value over polymorphic sites:\n")
print(merge(t1[, .(GROUP, Ho_published = round(Ho, 4), He_published = round(He, 4), Fis_published = round(Fis_vcftools, 4))],
            uni[SITES != "all called sites", .(GROUP, Ho_vcftools_style = round(HO_PER_ACCESSION_SAME_SITES, 4))], by = "GROUP"))
ua <- uni[SITES == "all called sites"]; up <- uni[SITES != "all called sites"]
ident <- intersect(c("N_oil", "oil_mean", "oil_sd", "n_oilseed", "n_confection", "n_restorer", "n_maintainer",
                     "n_breeding", "n_cultivar", "n_landrace", "top_region", "top_region_n"), names(t1))
t1new <- merge(t1[, c("GROUP", "N", ident), with = FALSE],
               ua[, .(GROUP, N_SNP_called = N_SITES, PctPoly = PCT_POLYMORPHIC, MAF, He = HE, Ho = HO, Fis = FIS, PIC)],
               by = "GROUP", all.x = TRUE)
t1new <- merge(t1new, up[, .(GROUP, N_polymorphic = N_SITES, He_polymorphic = HE, Ho_polymorphic = HO)], by = "GROUP", all.x = TRUE)
## the year each accession entered the collection, as a proxy for age
acc[, YEAR_NUM := suppressWarnings(as.integer(YEAR))]
yr <- acc[!is.na(YEAR_NUM), .(N_year = .N, year_q1 = as.integer(quantile(YEAR_NUM, .25)), year_median = as.integer(median(YEAR_NUM)),
                              year_q3 = as.integer(quantile(YEAR_NUM, .75))), by = GROUP]
t1new <- merge(t1new, yr, by = "GROUP", all.x = TRUE)
setcolorder(t1new, c("GROUP", "N", "N_SNP_called", "N_polymorphic", "PctPoly", "MAF", "He", "Ho", "Fis", "PIC",
                     "He_polymorphic", "Ho_polymorphic", ident, "N_year", "year_q1", "year_median", "year_q3"))
t1new <- t1new[order(match(GROUP, c(CL, "admixed")))]
print(rnd(copy(t1new)))
fwrite(t1new, file.path(TAB, "review_Table1_cluster_diversity_unified.csv"))

## ===========================================================================
## C. Heterozygous calls inside and outside the runs of homozygosity
##    Runs were called on imputed genotypes (step 07); the heterozygous calls
##    counted here are the non-imputed calls of the same accession.
## ===========================================================================
section("C. heterozygosity inside and outside runs of homozygosity")
fhom <- file.path(DV, "roh.hom")
if (file.exists(fhom)) {
  hom <- fread(fhom)
  sites <- data.table(CHROM = as.integer(pos$CHROM), START = pos$POS, END = pos$POS, IDX = seq_len(L))
  setkey(sites, CHROM, START, END)
  inout <- rbindlist(lapply(acc$SAMPLE, function(id) {
    seg <- hom[IID == id, .(CHROM = as.integer(CHR), START = as.integer(POS1), END = as.integer(POS2))]
    inr <- rep(FALSE, L)
    if (nrow(seg)) {
      ov <- foverlaps(seg, sites, type = "any", nomatch = NULL)
      inr[unique(ov$IDX)] <- TRUE
    }
    g <- G[id, ]; called <- !is.na(g); het <- called & g == 1L
    data.table(SAMPLE = id,
               CALLED_IN_ROH = sum(called & inr), CALLED_OUTSIDE_ROH = sum(called & !inr),
               HET_IN_ROH = sum(het & inr), HET_OUTSIDE_ROH = sum(het & !inr),
               HO_IN_ROH = sum(het & inr) / max(1, sum(called & inr)),
               HO_OUTSIDE_ROH = sum(het & !inr) / max(1, sum(called & !inr)),
               HO_ALL = sum(het) / sum(called),
               SHARE_OF_CALLED_SITES_IN_ROH = sum(called & inr) / sum(called),
               SHARE_OF_HET_CALLS_IN_ROH = sum(het & inr) / max(1, sum(het)))
  }))
  inout <- merge(inout, acc[, .(SAMPLE, GROUP, F_ROH, MISS_RATE)], by = "SAMPLE")
  fwrite(rnd(copy(inout)), file.path(TAB, "review_ho_roh_by_accession.csv"))
  byg <- inout[, .(N = .N, F_ROH = mean(F_ROH), SHARE_OF_CALLED_SITES_IN_ROH = mean(SHARE_OF_CALLED_SITES_IN_ROH),
                   HO_IN_ROH = mean(HO_IN_ROH), HO_OUTSIDE_ROH = mean(HO_OUTSIDE_ROH), HO_ALL = mean(HO_ALL),
                   RATIO_OUTSIDE_TO_INSIDE = mean(HO_OUTSIDE_ROH) / mean(HO_IN_ROH),
                   SHARE_OF_HET_CALLS_IN_ROH = mean(SHARE_OF_HET_CALLS_IN_ROH)),
               by = GROUP][order(match(GROUP, c(CL, "admixed")))]
  byg <- rbind(byg, inout[, .(GROUP = "All accessions", N = .N, F_ROH = mean(F_ROH),
                              SHARE_OF_CALLED_SITES_IN_ROH = mean(SHARE_OF_CALLED_SITES_IN_ROH),
                              HO_IN_ROH = mean(HO_IN_ROH), HO_OUTSIDE_ROH = mean(HO_OUTSIDE_ROH), HO_ALL = mean(HO_ALL),
                              RATIO_OUTSIDE_TO_INSIDE = mean(HO_OUTSIDE_ROH) / mean(HO_IN_ROH),
                              SHARE_OF_HET_CALLS_IN_ROH = mean(SHARE_OF_HET_CALLS_IN_ROH))])
  print(rnd(copy(byg)))
  fwrite(rnd(copy(byg)), file.path(TAB, "review_ho_roh_by_cluster.csv"))
} else message("no ", fhom)

## ===========================================================================
## D. Runs of homozygosity: without imputation, in the most complete libraries,
##    and by length class
## ===========================================================================
section("D. runs of homozygosity checks")
GENOME_KB <- sum(fread(file.path(ANA, "ref/chr_lengths.txt"), header = FALSE)$V2) / 1000
ri <- fread(file.path(DV, "roh.hom.indiv"))[, .(SAMPLE = IID, F_ROH_IMPUTED = KB / GENOME_KB, NSEG_IMPUTED = NSEG, KBAVG_IMPUTED = KBAVG)]
ri <- merge(ri, acc[, .(SAMPLE, GROUP, MISS_RATE)], by = "SAMPLE")
fnon <- file.path(RV, "roh_nonimputed.hom.indiv")
if (file.exists(fnon)) {
  rn <- fread(fnon)[, .(SAMPLE = IID, F_ROH_NONIMPUTED = KB / GENOME_KB, NSEG_NONIMPUTED = NSEG, KBAVG_NONIMPUTED = KBAVG)]
  rr <- merge(ri, rn, by = "SAMPLE")
  cat(sprintf("F_ROH imputed against non-imputed: Pearson r = %.3f, Spearman rho = %.3f\n",
              cor(rr$F_ROH_IMPUTED, rr$F_ROH_NONIMPUTED), cor(rr$F_ROH_IMPUTED, rr$F_ROH_NONIMPUTED, method = "spearman")))
  cat(sprintf("correlation with the missing rate: imputed r = %.3f, non-imputed r = %.3f\n",
              cor(rr$F_ROH_IMPUTED, rr$MISS_RATE), cor(rr$F_ROH_NONIMPUTED, rr$MISS_RATE)))
  byg <- rr[, .(N = .N, MISS_RATE = mean(MISS_RATE), F_ROH_IMPUTED = mean(F_ROH_IMPUTED), F_ROH_NONIMPUTED = mean(F_ROH_NONIMPUTED),
                NSEG_IMPUTED = mean(NSEG_IMPUTED), NSEG_NONIMPUTED = mean(NSEG_NONIMPUTED),
                KBAVG_IMPUTED = mean(KBAVG_IMPUTED), KBAVG_NONIMPUTED = mean(KBAVG_NONIMPUTED)),
            by = GROUP][order(-F_ROH_IMPUTED)]
  a1 <- anova(lm(F_ROH_NONIMPUTED ~ GROUP, data = rr))
  a2 <- anova(lm(F_ROH_NONIMPUTED ~ MISS_RATE + GROUP, data = rr))
  cat(sprintf("non-imputed F_ROH by cluster: F = %.2f on %d and %d df, p = %.2g; with the missing rate as covariate F = %.2f, p = %.2g\n",
              a1$`F value`[1], a1$Df[1], a1$Df[2], a1$`Pr(>F)`[1], a2$`F value`[2], a2$`Pr(>F)`[2]))
  print(rnd(copy(byg)))
  fwrite(rnd(copy(rr)), file.path(TAB, "review_roh_nonimputed_by_accession.csv"))
  fwrite(rnd(copy(byg)), file.path(TAB, "review_roh_nonimputed_by_cluster.csv"))
} else message("no ", fnon)

## the most complete libraries: missing rate at or below the panel median
med <- median(ri$MISS_RATE)
low <- ri[MISS_RATE <= med]
a3 <- anova(lm(F_ROH_IMPUTED ~ GROUP, data = low))
lowg <- low[, .(N = .N, MISS_RATE = mean(MISS_RATE), F_ROH = mean(F_ROH_IMPUTED)), by = GROUP][order(-F_ROH)]
lowg[, `:=`(SUBSET = sprintf("missing rate <= %.3f (panel median)", med), N_SUBSET = nrow(low),
            ANOVA_F = a3$`F value`[1], ANOVA_P = a3$`Pr(>F)`[1])]
allg <- ri[, .(GROUP, F_ROH_ALL_ACCESSIONS = F_ROH_IMPUTED)][, .(F_ROH_ALL_ACCESSIONS = mean(F_ROH_ALL_ACCESSIONS)), by = GROUP]
lowg <- merge(lowg, allg, by = "GROUP")[order(-F_ROH)]
print(rnd(copy(lowg)))
fwrite(rnd(copy(lowg)), file.path(TAB, "review_roh_low_missing_subset.csv"))

## length classes of the imputed-call runs
if (exists("hom")) {
  hom[, CLASS := cut(KB / 1000, c(0.5, 1, 2, 4, Inf), labels = c("0.5-1 Mb", "1-2 Mb", "2-4 Mb", ">4 Mb"), right = FALSE)]
  hom[, GROUP := acc$GROUP[match(IID, acc$SAMPLE)]]
  per_acc <- hom[, .(KB = sum(KB), NSEG = .N), by = .(IID, GROUP, CLASS)]
  grid <- CJ(IID = unique(hom$IID), CLASS = levels(hom$CLASS))
  grid[, GROUP := acc$GROUP[match(IID, acc$SAMPLE)]]
  per_acc <- merge(grid, per_acc, by = c("IID", "GROUP", "CLASS"), all.x = TRUE)
  per_acc[is.na(KB), `:=`(KB = 0, NSEG = 0L)]
  cls <- per_acc[, .(N = uniqueN(IID), MEAN_F_ROH_IN_CLASS = mean(KB) / GENOME_KB, MEAN_N_SEGMENTS = mean(NSEG)), by = .(GROUP, CLASS)]
  cls[, SHARE_OF_ROH_LENGTH := MEAN_F_ROH_IN_CLASS / sum(MEAN_F_ROH_IN_CLASS), by = GROUP]
  cls <- cls[order(match(GROUP, c(CL, "admixed")), CLASS)]
  print(rnd(copy(cls)))
  fwrite(rnd(copy(cls)), file.path(TAB, "review_roh_length_classes.csv"))
  cat("\nsegments by length class, whole panel:\n")
  print(hom[, .(SEGMENTS = .N, TOTAL_MB = sum(KB) / 1000, MEDIAN_LENGTH_KB = median(KB)), by = CLASS][order(CLASS)])
}

## ===========================================================================
## E. Clustering stability: ADMIXTURE against DAPC on all accessions,
##    agreement among the ten K = 5 runs, neighbouring K, and the run without
##    close relatives
## ===========================================================================
section("E. clustering stability")
ari <- function(a, b) {
  t <- table(a, b); n <- sum(t)
  comb2 <- function(x) x * (x - 1) / 2
  sij <- sum(comb2(t)); si <- sum(comb2(rowSums(t))); sj <- sum(comb2(colSums(t)))
  e <- si * sj / comb2(n)
  (sij - e) / ((si + sj) / 2 - e)
}
q  <- fread(find_tab("_admixture_Q\\.csv$"))
dp <- fread(find_tab("_dapc_coordinates\\.csv$"))
m  <- merge(q[, .(SAMPLE, CLUSTER, GROUP)], dp[, .(SAMPLE, DAPC)], by = "SAMPLE")
stab <- data.table(COMPARISON = c("ADMIXTURE cluster of highest ancestry against DAPC group, all accessions",
                                  "ADMIXTURE assignment (q >= 0.60) against DAPC group, assigned accessions",
                                  "ADMIXTURE assignment with admixed as a sixth label against DAPC group, all accessions"),
                   N = c(nrow(m), m[GROUP != "admixed", .N], nrow(m)),
                   ARI = c(ari(m$CLUSTER, m$DAPC), m[GROUP != "admixed", ari(GROUP, DAPC)], ari(m$GROUP, m$DAPC)))
cat("DAPC group by ADMIXTURE cluster of highest ancestry, all accessions:\n"); print(table(m$CLUSTER, m$DAPC))

fam  <- fread(file.path(PL, "sunflower_pruned.fam"), header = FALSE)$V2
best <- fread(file.path(AD, "best_rep_per_K.tsv"), header = FALSE, col.names = c("K", "REP", "CV"))
cv   <- fread(file.path(AD, "cv_summary.tsv"))
readQ <- function(K, rep) as.matrix(fread(file.path(AD, sprintf("K%d_rep%d/sunflower_pruned.%d.Q", K, rep, K))))
match_cols <- function(Qx, Qref) {
  cm <- cor(Qx, Qref); perm <- integer(ncol(Qref))
  for (i in seq_len(ncol(Qref))) {
    idx <- which(cm == max(cm), arr.ind = TRUE)[1, ]
    perm[idx[2]] <- idx[1]; cm[idx[1], ] <- -Inf; cm[, idx[2]] <- -Inf
  }
  Qx[, perm, drop = FALSE]
}
best5 <- best[K == 5, REP]
Qb <- readQ(5, best5)
stopifnot(nrow(Qb) == length(fam), identical(fam, q$SAMPLE[match(fam, q$SAMPLE)]))
base_cl  <- max.col(Qb)
base_grp <- ifelse(apply(Qb, 1, max) >= 0.60, base_cl, 0L)
cat(sprintf("\nreported K = 5 run: replicate %d (lowest cross-validation error, %.5f)\n", best5, best[K == 5, CV]))
reps <- rbindlist(lapply(cv[K == 5 & !is.na(CV_ERROR), REP], function(r) {
  Qr <- match_cols(readQ(5, r), Qb); cl <- max.col(Qr); grp <- ifelse(apply(Qr, 1, max) >= 0.60, cl, 0L)
  data.table(K = 5L, REP = r, CV_ERROR = cv[K == 5 & REP == r, CV_ERROR], REPORTED_RUN = r == best5,
             ARI_HIGHEST_ANCESTRY = ari(cl, base_cl), ARI_ASSIGNMENT = ari(grp, base_grp),
             ACCESSIONS_CHANGING_ASSIGNMENT = sum(grp != base_grp), MEAN_ABS_Q_DIFFERENCE = mean(abs(Qr - Qb)),
             MAX_ABS_Q_DIFFERENCE = max(abs(Qr - Qb)))
}))[order(REP)]
print(rnd(copy(reps), 5))
fwrite(rnd(copy(reps), 5), file.path(TAB, "review_admixture_k5_replicates.csv"))

nk <- rbindlist(lapply(c(4L, 6L), function(KK) {
  rep_k <- best[K == KK, REP][1]
  Qk <- readQ(KK, rep_k); cl <- max.col(Qk)
  cat(sprintf("\nK = %d (replicate %d) against K = 5, clusters of highest ancestry:\n", KK, rep_k))
  print(table(K5 = base_cl, other = cl))
  data.table(COMPARISON = sprintf("K = %d against K = 5, clusters of highest ancestry, all accessions", KK),
             N = length(cl), ARI = ari(cl, base_cl))
}))
stab <- rbind(stab, nk)

famn <- file.path(RV, "pruned_norel.fam")
if (file.exists(famn)) {
  keep <- fread(famn, header = FALSE)$V2; idx <- match(keep, fam)
  norel <- rbindlist(lapply(1:3, function(s) {
    Qn <- match_cols(as.matrix(fread(file.path(RV, sprintf("pruned_norel_K5_s%d.Q", s)))), Qb[idx, ])
    cl <- max.col(Qn); grp <- ifelse(apply(Qn, 1, max) >= 0.60, cl, 0L)
    data.table(SEED = s, N = length(keep),
               ARI_HIGHEST_ANCESTRY = ari(cl, base_cl[idx]), ARI_ASSIGNMENT = ari(grp, base_grp[idx]),
               ACCESSIONS_CHANGING_ASSIGNMENT = sum(grp != base_grp[idx]),
               ADMIXED_ORIGINAL = sum(base_grp[idx] == 0L), ADMIXED_WITHOUT_RELATIVES = sum(grp == 0L),
               MEAN_ABS_Q_DIFFERENCE = mean(abs(Qn - Qb[idx, ])),
               CLUSTER_SIZES = paste(tabulate(grp[grp > 0], 5), collapse = "/"))
  }))
  cat(sprintf("\nADMIXTURE at K = 5 without %d close relatives (%d accessions kept):\n", length(fam) - length(keep), length(keep)))
  print(rnd(copy(norel)))
  fwrite(rnd(copy(norel)), file.path(TAB, "review_admixture_without_relatives.csv"))
  stab <- rbind(stab, data.table(COMPARISON = "ADMIXTURE K = 5 without close relatives against the reported run, assignment (mean of three seeds)",
                                 N = length(keep), ARI = mean(norel$ARI_ASSIGNMENT)))
} else message("no ", famn)
print(rnd(copy(stab)))
fwrite(rnd(copy(stab)), file.path(TAB, "review_clustering_agreement.csv"))

## ===========================================================================
## F. Single-factor permutational analysis of variance on the identity-by-state
##    distance (plink --distance 1-ibs on the non-imputed set), 999 permutations
## ===========================================================================
section("F. single-factor PERMANOVA")
fmd <- file.path(RV, "ibs.mdist")
if (file.exists(fmd)) {
  D <- as.matrix(fread(fmd, header = FALSE)); dids <- fread(paste0(fmd, ".id"), header = FALSE)$V2
  dimnames(D) <- list(dids, dids)
  permanova <- function(ids, g, nperm = 999) {
    ok <- !is.na(g) & ids %in% dids; ids <- ids[ok]; g <- factor(g[ok])
    D2 <- D[ids, ids]^2; N <- length(g); a <- nlevels(g)
    sst <- sum(D2[upper.tri(D2)]) / N
    ssw <- function(g) { s <- 0
      for (l in levels(g)) { i <- which(g == l); if (length(i) > 1) s <- s + sum(D2[i, i][upper.tri(D2[i, i])]) / length(i) }
      s }
    w0 <- ssw(g); F0 <- ((sst - w0) / (a - 1)) / (w0 / (N - a))
    perm <- vapply(seq_len(nperm), function(k) { gp <- sample(g); wp <- ssw(gp); ((sst - wp) / (a - 1)) / (wp / (N - a)) }, numeric(1))
    data.table(N = N, LEVELS = a, PSEUDO_F = F0, R2 = (sst - w0) / sst, P = (sum(perm >= F0) + 1) / (nperm + 1))
  }
  reg5 <- acc[REGION != "Unknown", .N, by = REGION][N >= 10, REGION]
  set.seed(2026)
  pm <- rbindlist(list(
    cbind(FACTOR = "Region, all accessions in the five regions with at least ten", acc[REGION %in% reg5, permanova(SAMPLE, REGION)]),
    cbind(FACTOR = "Region, assigned accessions", acc[REGION %in% reg5 & GROUP != "admixed", permanova(SAMPLE, REGION)]),
    cbind(FACTOR = "Region, admixed accessions", acc[REGION %in% reg5 & GROUP == "admixed", permanova(SAMPLE, REGION)]),
    cbind(FACTOR = "Region, collected material", acc[REGION %in% reg5 & ORIGIN_EVENT == "Collected", permanova(SAMPLE, REGION)]),
    cbind(FACTOR = "Region within cultivars", acc[REGION %in% reg5 & STATUS == "Cultivar", permanova(SAMPLE, REGION)]),
    cbind(FACTOR = "Improvement status", acc[, permanova(SAMPLE, STATUS)]),
    cbind(FACTOR = "Market type (labeled accessions)", acc[!is.na(MARKET_TYPE), permanova(SAMPLE, MARKET_TYPE)]),
    cbind(FACTOR = "Heterotic group (labeled accessions)", acc[!is.na(HETEROTIC_GROUP), permanova(SAMPLE, HETEROTIC_GROUP)]),
    cbind(FACTOR = "Ancestry cluster, assigned accessions", acc[GROUP != "admixed", permanova(SAMPLE, GROUP)]),
    cbind(FACTOR = "Ancestry cluster with admixed as a sixth group", acc[, permanova(SAMPLE, GROUP)])))
  print(rnd(copy(pm)))
  fwrite(rnd(copy(pm)), file.path(TAB, "review_permanova.csv"))
} else message("no ", fmd)

## ===========================================================================
## G. The first principal component and seed oil content within groups
## ===========================================================================
section("G. PC1 and seed oil content")
pcs <- fread(find_tab("_pca_scores\\.csv$"))[, .(SAMPLE, PC1 = as.numeric(PC1))]
d <- merge(pcs, acc[, .(SAMPLE, OIL_PCT, STATUS, GROUP)], by = "SAMPLE")[!is.na(OIL_PCT)]
ct <- function(x, label) {
  if (nrow(x) < 8) return(NULL)
  t <- cor.test(x$PC1, x$OIL_PCT)
  data.table(SUBSET = label, N = nrow(x), R = unname(t$estimate), P = t$p.value)
}
res <- d[, lm(cbind(PC1, OIL_PCT) ~ GROUP)]$residuals
pooled <- cor.test(res[, 1], res[, 2])
pc <- rbindlist(c(list(ct(d, "All accessions with a recorded value"), ct(d[GROUP != "Q4"], "Without Q4"),
                       ct(d[GROUP != "admixed"], "Assigned accessions"), ct(d[GROUP == "admixed"], "Admixed accessions")),
                  lapply(c("Cultivar", "Breeding line", "Landrace"), function(s) ct(d[STATUS == s], s)),
                  lapply(CL, function(g) ct(d[GROUP == g], paste("Within", g))),
                  list(data.table(SUBSET = "Pooled within the six groups (residuals from group means)", N = nrow(d),
                                  R = unname(pooled$estimate), P = pooled$p.value))))
m0 <- lm(OIL_PCT ~ GROUP, d); m1 <- lm(OIL_PCT ~ PC1, d); m2 <- lm(OIL_PCT ~ GROUP + PC1, d)
an <- anova(m0, m2)
cat(sprintf("R2 of seed oil content: PC1 alone %.3f, group alone %.3f, group and PC1 %.3f; PC1 added to group: F = %.1f, p = %.2g\n",
            summary(m1)$r.squared, summary(m0)$r.squared, summary(m2)$r.squared, an$F[2], an$`Pr(>F)`[2]))
pc <- rbind(pc, data.table(SUBSET = c("R2, PC1 alone", "R2, group alone", "R2, group and PC1", "PC1 added to group (partial F test)"),
                           N = nrow(d), R = c(summary(m1)$r.squared, summary(m0)$r.squared, summary(m2)$r.squared, an$F[2]),
                           P = c(NA, NA, NA, an$`Pr(>F)`[2])))
print(rnd(copy(pc), 5))
fwrite(rnd(copy(pc), 5), file.path(TAB, "review_pc1_oil.csv"))

## ===========================================================================
## H. Multiple-testing correction of the Kolmogorov-Smirnov tests on the cores
## ===========================================================================
section("H. Kolmogorov-Smirnov tests corrected for the number of descriptors")
cov <- fread(find_tab("_core_phenotype_coverage\\.csv$"))
ks <- cov[TYPE == "continuous" & !is.na(KS_P), .(DESCRIPTOR, DESCRIPTION, SCHEME, CORE_N, PANEL_N, KS_P)]
ks[, `:=`(N_TESTS = .N, KS_P_HOLM = p.adjust(KS_P, "holm"), KS_P_BH = p.adjust(KS_P, "BH")), by = SCHEME]
print(ks)
fwrite(ks, file.path(TAB, "review_core_ks_adjusted.csv"))

## ===========================================================================
## I. Share of the assembled chromosomes within a given distance of a marker
## ===========================================================================
section("I. genome coverage of the marker set")
bim <- fread(file.path(PL, "sunflower.bim"), header = FALSE)[, .(CHROM = as.integer(V1), POS = as.numeric(V4))]
len <- fread(file.path(ANA, "ref/chr_lengths.txt"), header = FALSE)[, .(CHROM = as.integer(V1), LEN = as.numeric(V2))]
covered <- function(p, Lc, d) {
  s <- pmax(1, p - d); e <- pmin(Lc, p + d); o <- order(s); s <- s[o]; e <- e[o]
  cm <- cummax(e); new <- c(TRUE, s[-1] > cm[-length(cm)])
  starts <- s[new]; ends <- c(cm[which(new)[-1] - 1], cm[length(cm)])
  sum(ends - starts + 1)
}
gc <- rbindlist(lapply(c(5, 10, 25, 50, 100, 200, 500) * 1000, function(d) {
  cv <- sum(vapply(len$CHROM, function(c) covered(bim[CHROM == c, POS], len[CHROM == c, LEN], d), numeric(1)))
  data.table(DISTANCE_KB = d / 1000, SHARE_OF_GENOME_WITHIN_DISTANCE = cv / sum(len$LEN))
}))
gaps <- bim[order(CHROM, POS), .(GAP = diff(POS)), by = CHROM]$GAP
cat(sprintf("marker spacing: median %.2f kb, mean %.1f kb, 90th percentile %.1f kb, 99th percentile %.0f kb, largest gap %.2f Mb\n",
            median(gaps) / 1000, mean(gaps) / 1000, quantile(gaps, .9) / 1000, quantile(gaps, .99) / 1000, max(gaps) / 1e6))
gc <- rbind(gc, data.table(DISTANCE_KB = c(NA, NA, NA, NA),
                           SHARE_OF_GENOME_WITHIN_DISTANCE = c(median(gaps), mean(gaps), quantile(gaps, .9), max(gaps)) / 1000))
gc[, NOTE := c(rep("share of assembled chromosome length within the distance of a retained marker", 7),
               "median spacing between adjacent markers, kb", "mean spacing, kb", "90th percentile of spacing, kb", "largest gap, kb")]
print(rnd(copy(gc)))
fwrite(rnd(copy(gc)), file.path(TAB, "review_genome_coverage.csv"))

## ===========================================================================
## J. Linkage disequilibrium at a matched sample size: ten random subsets of
##    23 accessions, summarized as in R/15
## ===========================================================================
section("J. linkage disequilibrium in random subsets of 23 accessions")
read_ld <- function(f) {
  d <- if (grepl("\\.gz$", f)) fread(cmd = paste("zcat", shQuote(f))) else fread(f)
  setnames(d, 1:2, c("Dist", "MeanR2"))
  np <- if (ncol(d) >= 6) d[[6]] else rep(1L, nrow(d))
  data.table(dist_bp = as.numeric(d$Dist), r2 = as.numeric(d$MeanR2), n = as.numeric(np))
}
BIN <- 2000
bin_curve <- function(x) {
  e <- x[, .(r2 = sum(r2 * n) / sum(n), pairs = sum(n)), by = .(SET, bin = floor(dist_bp / BIN))]
  e[, dist_kb := (bin + 0.5) * BIN / 1000]; setorder(e, SET, bin); e[]
}
fixed_points <- function(e) e[, {
  o <- order(bin); ir <- isoreg(bin[o], -r2[o]); mono <- -ir$yf; d <- (bin[o] + 0.5) * BIN / 1000
  at <- function(x) mono[which.min(abs(d - x))]
  f <- function(t) { i <- which(mono <= t)[1]; if (is.na(i)) NA_real_ else d[i] }
  .(R2_1KB = at(1), R2_5KB = at(5), R2_10KB = at(10), R2_20KB = at(20), R2_50KB = at(50), R2_100KB = at(100),
    R2_200KB = at(200), R2_500KB = at(500), R2_800KB_TO_1MB = mean(r2[o][d >= 800]),
    DIST_R2_0.20_KB = f(0.20), DIST_R2_0.15_KB = f(0.15), DIST_R2_0.10_KB = f(0.10))
}, by = SET]
subs <- list.files(RV, pattern = "^sub23_[0-9]+\\.stat(\\.gz)?$", full.names = TRUE)
if (length(subs)) {
  lds <- rbindlist(lapply(subs, function(f) { d <- read_ld(f); d[, SET := sub("\\.stat(\\.gz)?$", "", basename(f))]; d }))
  lds <- lds[dist_bp > 0 & is.finite(r2)]
  fp <- fixed_points(bin_curve(lds))
  cols <- setdiff(names(fp), "SET")
  summ <- rbind(fp[, c(SET = "Mean of ten subsets of 23", lapply(.SD, mean, na.rm = TRUE)), .SDcols = cols],
                fp[, c(SET = "SD across subsets", lapply(.SD, sd, na.rm = TRUE)), .SDcols = cols],
                fp[, c(SET = "Minimum", lapply(.SD, min, na.rm = TRUE)), .SDcols = cols],
                fp[, c(SET = "Maximum", lapply(.SD, max, na.rm = TRUE)), .SDcols = cols])
  ## the published curves, recomputed here with the same code so the background
  ## level between 800 kb and 1 Mb is available for every set
  main <- list.files(LD, pattern = "^(all|grp_Q[0-9]+)\\.stat(\\.gz)?$", full.names = TRUE)
  ldm <- rbindlist(lapply(main, function(f) { d <- read_ld(f)
    d[, SET := sub("^grp_", "", sub("^all$", "All accessions", sub("\\.stat(\\.gz)?$", "", basename(f))))]; d }))
  fpm <- fixed_points(bin_curve(ldm[dist_bp > 0 & is.finite(r2)]))
  out <- rbind(fpm, summ, fp, fill = TRUE)
  print(rnd(copy(out), 3))
  fwrite(rnd(copy(out), 4), file.path(TAB, "review_ld_subsets_of_23.csv"))
  ## binned curves of the subsets, for the ratio panel of the LD figure (R/20)
  eb <- bin_curve(lds)
  fwrite(eb[, .(SET, DISTANCE_KB = dist_kb, R2 = round(r2, 5), PAIRS = pairs)],
         file.path(TAB, "review_ld_subsets_binned.csv"))
} else message("no sub23 stat files in ", RV)

## ===========================================================================
## K. Windowed F_ST between Q2 and Q5 (step 05, 100-kb windows)
## ===========================================================================
section("K. windowed F_ST scan")
wf <- list.files(DV, pattern = "^fst_Q[0-9]_vs_Q[0-9]\\.windowed\\.weir\\.fst$", full.names = TRUE)
if (length(wf)) {
  fscan <- rbindlist(lapply(wf, function(f) {
    w <- fread(f)[N_VARIANTS >= 10]
    pair <- sub("^fst_(Q[0-9])_vs_(Q[0-9]).*$", "\\1-\\2", basename(f))
    thr <- quantile(w$WEIGHTED_FST, 0.99)
    top <- w[WEIGHTED_FST >= thr]
    data.table(PAIR = pair, WINDOWS = nrow(w), MEDIAN_FST = median(w$WEIGHTED_FST), P99_FST = thr, MAX_FST = max(w$WEIGHTED_FST),
               TOP_1PCT_WINDOWS = nrow(top), CHROMOSOMES_WITH_TOP_WINDOWS = uniqueN(top$CHROM),
               TOP_WINDOWS_BY_CHROMOSOME = paste(sprintf("Ha%s:%d", names(table(top$CHROM)), as.integer(table(top$CHROM))), collapse = " "))
  }))
  print(rnd(copy(fscan), 3))
  fwrite(rnd(copy(fscan), 4), file.path(TAB, "review_fst_scan_summary.csv"))
  f25 <- file.path(DV, "fst_Q2_vs_Q5.windowed.weir.fst")
  if (file.exists(f25)) {
    w <- fread(f25)[N_VARIANTS >= 10]
    w[, TOP_1PCT := WEIGHTED_FST >= quantile(WEIGHTED_FST, 0.99)]
    setorder(w, -WEIGHTED_FST)
    cat("\ntwenty most differentiated 100-kb windows, Q2 against Q5:\n")
    print(w[1:20, .(CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST = round(WEIGHTED_FST, 3))])
    fwrite(w, file.path(TAB, "review_fst_windows_Q2_Q5.csv"))
  }
} else message("no windowed F_ST files in ", DV)

## ===========================================================================
## L. KING kinship of the three libraries of PI 490282 (276-library set)
## ===========================================================================
section("L. KING kinship of the technical replicates")
fk <- file.path(RV, "king276.kin0")
if (file.exists(fk)) {
  k <- fread(fk)
  setnames(k, sub("^#", "", names(k)))
  # plink2 writes IID1/IID2 (or FID1/IID1/FID2/IID2 when family ids are present)
  if ("IID1" %in% names(k)) setnames(k, c("IID1", "IID2"), c("ID1", "ID2"))
  # both members must be libraries of PI 490282; other accessions that pass the
  # table filter with one of them are not replicates
  rep <- k[grepl("^PI_490282", ID1) & grepl("^PI_490282", ID2)]
  dup <- k[(ID1 %in% c("PI_599764", "PI_497931") & ID2 %in% c("PI_599764", "PI_497931")) |
           (ID1 %in% c("PI_664196", "PI_664191") & ID2 %in% c("PI_664196", "PI_664191"))]
  out <- rbind(cbind(SET = "replicate libraries of PI 490282", rep), cbind(SET = "candidate duplicate accessions", dup))
  print(out)
  fwrite(out, file.path(TAB, "review_king_replicates.csv"))
  cat(sprintf("pairs with kinship >= 0.354 among the 276 libraries: %d\n", k[KINSHIP >= 0.354, .N]))
} else message("no ", fk)

## ===========================================================================
## M. F_ST, d_xy and Nei's distance for predefined groupings on the same
##    accessions, from the per-site Weir and Cockerham components
## ===========================================================================
section("M. differentiation of predefined groupings and absolute divergence")
wc_pair <- function(ids1, ids2) {
  s1 <- site_stats(ids1); s2 <- site_stats(ids2)
  n1 <- s1$n; n2 <- s2$n; ok <- n1 >= 2 & n2 >= 2
  r <- 2; nbar <- (n1 + n2) / 2; nsum <- n1 + n2
  nc <- (nsum - (n1^2 + n2^2) / nsum) / (r - 1)
  pbar <- (n1 * s1$p + n2 * s2$p) / nsum
  s2v <- (n1 * (s1$p - pbar)^2 + n2 * (s2$p - pbar)^2) / ((r - 1) * nbar)
  hbar <- (n1 * s1$ho + n2 * s2$ho) / nsum
  a <- nbar / nc * (s2v - (pbar * (1 - pbar) - (r - 1) / r * s2v - hbar / 4) / (nbar - 1))
  b <- nbar / (nbar - 1) * (pbar * (1 - pbar) - (r - 1) / r * s2v - (2 * nbar - 1) / (4 * nbar) * hbar)
  c <- hbar / 2
  p1 <- s1$p; p2 <- s2$p
  list(ok = ok, a = a, abc = a + b + c,
       dxy = p1 * (1 - p2) + p2 * (1 - p1),
       jx = p1^2 + (1 - p1)^2, jy = p2^2 + (1 - p2)^2, jxy = p1 * p2 + (1 - p1) * (1 - p2))
}
pair_row <- function(set, g1, g2, ids1, ids2) {
  w <- wc_pair(ids1, ids2); o <- w$ok
  data.table(SET = set, GROUP1 = g1, GROUP2 = g2, N1 = length(ids1), N2 = length(ids2), N_SITES = sum(o),
             FST = sum(w$a[o]) / sum(w$abc[o]), DXY = mean(w$dxy[o]),
             NEI_D = -log(mean(w$jxy[o]) / sqrt(mean(w$jx[o]) * mean(w$jy[o]))))
}
all_pairs <- function(set, dt, col, min_n) {
  lv <- dt[, .N, by = c(col)][N >= min_n][[col]]
  if (length(lv) < 2) return(NULL)
  cmb <- combn(sort(lv), 2)
  rbindlist(lapply(seq_len(ncol(cmb)), function(i)
    pair_row(set, cmb[1, i], cmb[2, i], dt[get(col) == cmb[1, i], SAMPLE], dt[get(col) == cmb[2, i], SAMPLE])))
}
reg <- acc[REGION != "Unknown"]
fst <- rbindlist(list(
  all_pairs("regions, all accessions", reg, "REGION", 10),
  all_pairs("regions, assigned accessions", reg[GROUP != "admixed"], "REGION", 10),
  all_pairs("regions, admixed accessions", reg[GROUP == "admixed"], "REGION", 10),
  all_pairs("regions, collected material (at least five accessions)", reg[ORIGIN_EVENT == "Collected"], "REGION", 5),
  all_pairs("regions within cultivars", reg[STATUS == "Cultivar"], "REGION", 10),
  all_pairs("regions within breeding lines", reg[STATUS == "Breeding line"], "REGION", 10),
  rbindlist(lapply(CL, function(g) all_pairs(paste("regions within cluster", g), reg[GROUP == g], "REGION", 10))),
  all_pairs("improvement status, all accessions", acc, "STATUS", 10),
  all_pairs("heterotic group, labeled accessions", acc[!is.na(HETEROTIC_GROUP)], "HETEROTIC_GROUP", 5),
  all_pairs("market type, labeled accessions", acc[!is.na(MARKET_TYPE)], "MARKET_TYPE", 5),
  all_pairs("ancestry clusters, assigned accessions", acc[GROUP != "admixed"], "GROUP", 5),
  rbindlist(lapply(CL, function(g) pair_row("ancestry cluster against the admixed accessions", g, "admixed",
                                           acc[GROUP == g, SAMPLE], acc[GROUP == "admixed", SAMPLE])))))
print(rnd(copy(fst)))
fwrite(rnd(copy(fst), 5), file.path(TAB, "review_fst_predefined_groups.csv"))
cat("\nsummary by set:\n")
print(fst[, .(PAIRS = .N, MIN_FST = round(min(FST), 4), MEDIAN_FST = round(median(FST), 4), MAX_FST = round(max(FST), 4)), by = SET])
cat("\nvalidation against the published pairwise F_ST (Table S10):\n")
t10 <- fread(find_tab("_fst_intervals\\.csv$"))
print(merge(t10[, .(GROUP1 = POP1, GROUP2 = POP2, PUBLISHED = FST)],
            fst[SET == "ancestry clusters, assigned accessions", .(GROUP1, GROUP2, RECOMPUTED = round(FST, 4))], by = c("GROUP1", "GROUP2")))
cl_div <- fst[SET == "ancestry clusters, assigned accessions"][order(-FST)]
cl_div[, `:=`(RANK_FST = rank(-FST), RANK_DXY = rank(-DXY), RANK_NEI_D = rank(-NEI_D))]
print(rnd(copy(cl_div)))
fwrite(rnd(copy(cl_div), 5), file.path(TAB, "review_divergence_clusters.csv"))

## Q2-Q5 against Q2-Q4: is the difference larger than its uncertainty?
w25 <- wc_pair(acc[GROUP == "Q2", SAMPLE], acc[GROUP == "Q5", SAMPLE])
w24 <- wc_pair(acc[GROUP == "Q2", SAMPLE], acc[GROUP == "Q4", SAMPLE])
chrom <- as.integer(pos$CHROM)
ratio <- function(w, sel) sum(w$a[sel & w$ok]) / sum(w$abc[sel & w$ok])
full <- ratio(w25, TRUE) - ratio(w24, TRUE)
chrs <- sort(unique(chrom))
loo <- vapply(chrs, function(c) ratio(w25, chrom != c) - ratio(w24, chrom != c), numeric(1))
k <- length(chrs); se <- sqrt((k - 1) / k * sum((loo - mean(loo))^2))
per_chr <- vapply(chrs, function(c) ratio(w25, chrom == c) - ratio(w24, chrom == c), numeric(1))
set.seed(7)
boot <- replicate(1000, { i <- sample.int(L, L, replace = TRUE); sel <- tabulate(i, L)
  (sum((w25$a * sel)[w25$ok]) / sum((w25$abc * sel)[w25$ok])) - (sum((w24$a * sel)[w24$ok]) / sum((w24$abc * sel)[w24$ok])) })
diff <- data.table(STATISTIC = c("F_ST Q2-Q5", "F_ST Q2-Q4", "difference", "jackknife SE of the difference (17 chromosomes)",
                                 "z of the difference", "chromosomes where Q2-Q5 exceeds Q2-Q4", "bootstrap over loci, 2.5th percentile",
                                 "bootstrap over loci, 97.5th percentile"),
                   VALUE = c(ratio(w25, TRUE), ratio(w24, TRUE), full, se, full / se, sum(per_chr > 0), quantile(boot, .025), quantile(boot, .975)))
print(rnd(copy(diff)))
fwrite(rnd(copy(diff), 5), file.path(TAB, "review_fst_q2q5_vs_q2q4.csv"))

## ===========================================================================
## O. Membership of the two 20 percent cores, and the duplicate pairs in them
## ===========================================================================
section("O. core membership")
CR <- file.path(OUT, "07_core")
f_c <- file.path(CR, c("core_20.txt", "core_stratified_20.txt"))
if (all(file.exists(f_c))) {
  cw <- readLines(f_c[1]); cs <- readLines(f_c[2])
  mem <- acc[, .(SAMPLE, GROUP, IN_WHOLE_PANEL_CORE = as.integer(SAMPLE %in% cw), IN_STRATIFIED_CORE = as.integer(SAMPLE %in% cs))]
  cat(sprintf("whole-panel core %d accessions, stratified core %d\n", sum(mem$IN_WHOLE_PANEL_CORE), sum(mem$IN_STRATIFIED_CORE)))
  for (pr in list(c("PI_599764", "PI_497931"), c("PI_664196", "PI_664191")))
    cat(sprintf("duplicate pair %s / %s: whole-panel core %d/%d, stratified core %d/%d\n", pr[1], pr[2],
                pr[1] %in% cw, pr[2] %in% cw, pr[1] %in% cs, pr[2] %in% cs))
  fwrite(mem, file.path(TAB, "review_core_membership.csv"))
} else message("no core lists in ", CR)

## ===========================================================================
## N. The sites removed by the heterozygous call rate filter, within clusters
##    (needs results/12_review/allsites.012 from slurm/14)
## ===========================================================================
section("N. removed sites within clusters")
fall <- file.path(RV, "allsites.012")
if (file.exists(fall)) {
  rm(G); gc()
  aind <- readLines(paste0(fall, ".indv"))
  apos <- fread(paste0(fall, ".pos"), header = FALSE, col.names = c("CHROM", "POS"))
  GA <- matrix(scan(fall, what = integer(), sep = "\t", quiet = TRUE), nrow = length(aind), byrow = TRUE)[, -1, drop = FALSE]
  GA[GA == -1L] <- NA_integer_; rownames(GA) <- aind
  GA <- GA[acc$SAMPLE, , drop = FALSE]
  n <- colSums(!is.na(GA)); het <- colSums(GA == 1L, na.rm = TRUE) / pmax(n, 1)
  p <- colSums(GA, na.rm = TRUE) / pmax(2 * n, 1); maf <- pmin(p, 1 - p)
  base <- n / nrow(GA) >= 0.5 & maf >= 0.05
  removed <- base & het > 0.55; kept <- base & het <= 0.55
  cat(sprintf("%d sites in the export; %d pass call rate and MAF; %d removed by the heterozygous call rate ceiling; %d retained\n",
              ncol(GA), sum(base), sum(removed), sum(kept)))
  rs <- rbindlist(lapply(c(CL, "admixed"), function(g) {
    ids <- acc[GROUP == g, SAMPLE]
    rate <- function(mask) { sub <- GA[ids, mask, drop = FALSE]; nn <- colSums(!is.na(sub)); h <- colSums(sub == 1L, na.rm = TRUE) / pmax(nn, 1); h[nn > 0] }
    hr <- rate(removed); hk <- rate(kept)
    data.table(GROUP = g, N = length(ids), REMOVED_SITES_CALLED = length(hr), HET_RATE_REMOVED_SITES = mean(hr),
               SHARE_REMOVED_SITES_HET_ABOVE_0_5 = mean(hr > 0.5), HET_RATE_RETAINED_SITES = mean(hk))
  }))
  print(rnd(copy(rs)))
  fwrite(rnd(copy(rs)), file.path(TAB, "review_removed_sites_by_cluster.csv"))
} else message("no ", fall, " (run slurm/14_review_extra.slurm)")

cat("step 24 done\n")
