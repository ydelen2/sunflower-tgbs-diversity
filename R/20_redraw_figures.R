## Step 20 - manuscript figures drawn from saved results ---------------------------
## Every panel is drawn from files written by steps 10 to 16; no analysis is
## recomputed.
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))

FIG <- file.path(OUT, "08_figures_redraw"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
QC  <- file.path(OUT, "01_qc")
AD  <- file.path(OUT, "03_admixture")

have <- function(f) {
  ok <- file.exists(f)
  if (!ok) message("not found, skipped: ", f)
  ok
}

## a results table, found by the end of its file name
find_tab <- function(pattern, dir = TAB) {
  f <- list.files(dir, pattern, full.names = TRUE)
  if (length(f)) f[1] else file.path(dir, pattern)
}

## Base font size of each figure. Figures are reduced by different factors when
## they are placed on the page, so the base size is set per figure such that
## axis titles print at about 8 pt and tick labels at about 6.5 pt.
FS <- c(fig1ab = 15, fig1c = 14, fig2 = 13, fig3 = 14.5, fig4 = 11, fig5 = 13.5,
        s1 = 7.5, s2 = 9, s3 = 11, s4 = 15, s6 = 7.5, s7 = 11, s8 = 13.5,
        s9 = 8.5, s10 = 13.5, s11 = 12.5)

## cluster order and key labels shared by every figure
GRP_LEVELS <- c("Q1", "Q2", "Q3", "Q4", "Q5", "admixed")
grp_label  <- function(x) sub("^admixed$", "Admixed", x)
scale_cluster <- function(aes_name = "colour", ...)
  scale_discrete_manual(aes_name, values = CLUSTER_COLS, breaks = GRP_LEVELS,
                        labels = grp_label, ...)

## ---- Figure S4: quality control histograms ------------------------------------
f_ss <- file.path(QC, "sample_stats.tsv")
if (have(f_ss)) {
  ss_all <- fread(f_ss)
  thr <- as.numeric(Sys.getenv("MAX_SAMPLE_MISS", "0.90"))
  p <- ggplot(ss_all, aes(MISS_RATE)) +
    geom_histogram(bins = 40, fill = PAL[1], colour = "white") +
    geom_vline(xintercept = thr, linetype = 2, colour = PAL[4]) +
    labs(x = "Missing genotype rate per library", y = "Number of libraries") + theme_pub(FS[["s4"]])
  save_fig(p, "FigS2_sample_missingness", 6, 4)
}

f_st <- file.path(QC, "site_stats_filtered.tsv")
if (have(f_st)) {
  st <- fread(f_st)[PASS == 1]
  st[, CHRL := factor(CHR_LAB[as.character(CHROM)], levels = CHR_LAB)]
  p <- ggplot(st, aes(MAF)) + geom_histogram(binwidth = .01, fill = PAL[3], colour = "white") +
    labs(x = "Minor allele frequency", y = "Number of SNPs") + theme_pub(FS[["s4"]])
  save_fig(p, "FigS3_maf_spectrum", 6, 4)
  p <- ggplot(st, aes(MISS)) + geom_histogram(binwidth = .02, fill = PAL[2], colour = "white") +
    labs(x = "Missing rate per SNP", y = "Number of SNPs") + theme_pub(FS[["s4"]])
  save_fig(p, "FigS4_site_missing", 6, 4)
  p <- ggplot(st, aes(HO)) + geom_histogram(binwidth = .01, fill = PAL[5], colour = "white") +
    labs(x = "Heterozygous call rate per SNP", y = "Number of SNPs") + theme_pub(FS[["s4"]])
  save_fig(p, "FigS5_site_het", 6, 4)

  ## Figure S3: each 1-Mb window drawn as a rectangle from its start to its end.
  ## The small overlap closes the hairline gaps that anti-aliasing leaves
  ## between adjacent tiles.
  dens <- st[!is.na(CHRL), .N, by = .(CHRL, WIN = floor(POS / 1e6))]
  dens[, Y := 18 - as.integer(CHRL)]
  p <- ggplot(dens) +
    geom_rect(aes(xmin = WIN, xmax = WIN + 1.08, ymin = Y - .4, ymax = Y + .4, fill = N)) +
    scale_fill_viridis_c(name = "SNPs\nper Mb", option = "C") +
    scale_y_continuous(breaks = 17:1, labels = CHR_LAB, expand = c(0, .3)) +
    scale_x_continuous(expand = c(.01, 0)) +
    labs(x = "Position (Mb)", y = NULL) + theme_pub(FS[["s3"]]) +
    theme(panel.grid.major = element_blank())
  save_fig(p, "Fig1_snp_density", 9, 5)
}

## ---- Figure S11: heterozygosity against missing rate --------------------------
f_sk <- file.path(QC, "sample_stats_kept.tsv")
if (have(f_sk)) {
  ss <- fread(f_sk)
  ss[, HET_PER_ALL := HO * (1 - MISS_RATE)]
  pl <- melt(ss[, .(SAMPLE, MISS_RATE,
                    `(a) Heterozygous / called sites` = HO,
                    `(b) Heterozygous / called, common sites` = CORE_HO,
                    `(c) Heterozygous / all sites` = HET_PER_ALL)],
             id.vars = c("SAMPLE", "MISS_RATE"), variable.name = "which", value.name = "Ho")
  p <- ggplot(pl, aes(MISS_RATE, Ho)) +
    geom_point(alpha = .6, size = 1.6, colour = PAL[1]) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = PAL[4], linewidth = .6) +
    facet_wrap(~which) +
    labs(x = "Missing genotype rate per accession", y = "Observed heterozygosity") +
    theme_pub(FS[["s11"]])
  save_fig(p, "FigS1_het_vs_missing", 10, 4)
}

## ---- Figure 1a and 1c, and Figure S8: ADMIXTURE -------------------------------
f_cv <- file.path(AD, "cv_summary.tsv")
if (have(f_cv)) {
  cv  <- fread(f_cv)[!is.na(CV_ERROR)]
  agg <- cv[, .(mean = mean(CV_ERROR), sd = sd(CV_ERROR)), by = K][order(K)]
  bestK <- agg[which.min(mean), K]
  p <- ggplot(agg, aes(K, mean)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .2, colour = "grey60") +
    geom_line(colour = PAL[1]) + geom_point(size = 2.4, colour = PAL[1]) +
    geom_point(data = agg[K == bestK], colour = PAL[4], size = 3.6) +
    scale_x_continuous(breaks = agg$K) +
    labs(x = "Number of ancestral clusters (K)",
         y = sprintf("Cross-validation error\n(mean ± SD, %d runs)", max(cv$REP))) +
    theme_pub(FS[["fig1ab"]])
  save_fig(p, "Fig3a_admixture_cv", 6.5, 4)
}

f_best <- file.path(AD, "best_rep_per_K.tsv")
f_fam  <- file.path(QC, "sunflower_pruned.fam")
if (have(f_best) && have(f_fam) && exists("bestK")) {
  ids  <- fread(f_fam, header = FALSE)$V2
  best <- fread(f_best, header = FALSE, col.names = c("K", "REP", "CV"))[K >= 2][order(K)]
  qlist <- lapply(seq_len(nrow(best)), function(i) {
    K <- best$K[i]
    f <- file.path(AD, sprintf("K%d_rep%d", K, best$REP[i]), sprintf("sunflower_pruned.%d.Q", K))
    q <- as.matrix(fread(f)); rownames(q) <- ids; q
  })
  names(qlist) <- best$K

  ## Match the columns of neighboring K by the correlation of their ancestry
  ## proportions, starting from the best K, so a cluster keeps its color from
  ## row to row. Above the best K each new cluster takes the next free color.
  match_cols <- function(A, B) {
    C <- cor(A, B); m <- rep(NA_integer_, ncol(A))
    for (s in seq_len(min(dim(C)))) {
      w <- which(C == max(C, na.rm = TRUE), arr.ind = TRUE)[1, ]
      m[w[1]] <- w[2]; C[w[1], ] <- NA; C[, w[2]] <- NA
    }
    m
  }
  Ks <- as.integer(names(qlist))
  cols <- list(); cols[[as.character(bestK)]] <- PAL[seq_len(bestK)]
  for (K in rev(Ks[Ks < bestK]))
    cols[[as.character(K)]] <- cols[[as.character(K + 1)]][
      match_cols(qlist[[as.character(K)]], qlist[[as.character(K + 1)]])]
  used <- bestK
  for (K in Ks[Ks > bestK]) {
    cc <- cols[[as.character(K - 1)]][
      match_cols(qlist[[as.character(K)]], qlist[[as.character(K - 1)]])]
    for (j in which(is.na(cc))) { used <- used + 1; cc[j] <- PAL[used] }
    cols[[as.character(K)]] <- cc
  }

  Qb  <- qlist[[as.character(bestK)]]
  ord <- order(paste0("Q", max.col(Qb)), -apply(Qb, 1, max))
  long <- rbindlist(lapply(names(qlist), function(k) {
    q <- qlist[[k]]
    data.table(SAMPLE = rep(ids, ncol(q)), K = as.integer(k),
               FILL = rep(cols[[k]], each = length(ids)), Q = as.vector(q))
  }))
  long[, SAMPLE := factor(SAMPLE, levels = ids[ord])]

  p <- ggplot(long, aes(SAMPLE, Q, fill = FILL)) +
    geom_col(width = 1) +
    facet_grid(K ~ ., switch = "y", labeller = as_labeller(function(x) paste("K =", x))) +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    labs(x = sprintf("%d accessions", length(ids)), y = "Ancestry proportion") +
    theme_pub(FS[["s8"]]) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.spacing.y = unit(12, "pt"),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
  save_fig(p, "Fig3b_admixture_barplots_K2toKmax", 11, 1 + 0.9 * length(qlist))

  kb <- long[K == bestK]
  kb[, CL := factor(FILL, levels = PAL[seq_len(bestK)], labels = paste0("Q", seq_len(bestK)))]
  p <- ggplot(kb, aes(SAMPLE, Q, fill = CL)) +
    geom_col(width = 1) +
    scale_fill_manual(values = setNames(PAL[seq_len(bestK)], paste0("Q", seq_len(bestK))),
                      name = "Cluster") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Ancestry proportion", title = sprintf("K = %d", bestK)) +
    theme_pub(FS[["fig1c"]]) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  save_fig(p, sprintf("Fig3c_admixture_K%d", bestK), 11, 3)
}

## ---- Figure 1b: discriminant analysis, at the height of panel a ---------------
f_dp <- find_tab("dapc_coordinates\\.csv$")
if (have(f_dp)) {
  dt <- fread(f_dp)
  p <- ggplot(dt, aes(LD1, LD2)) +
    stat_ellipse(aes(group = DAPC), level = .67, linewidth = .4, colour = "grey40") +
    geom_point(aes(colour = GROUP), size = 2.2, alpha = .85) +
    scale_cluster(name = "Cluster") +
    labs(x = "Discriminant function 1", y = "Discriminant function 2") + theme_pub(FS[["fig1ab"]])
  save_fig(p, "Fig5_dapc", 5.6, 4)
}

## ---- Figure 2, Figures S1, S6 and S9: principal components and oil content ----
f_pc <- find_tab("pca_scores\\.csv$")
f_pv <- find_tab("pca_variance\\.csv$")
if (have(f_pc) && have(f_pv)) {
  vec <- fread(f_pc)
  pct <- fread(f_pv)$VAR_PCT
  xl <- sprintf("PC1 (%.2f%%)", pct[1]); yl <- sprintf("PC2 (%.2f%%)", pct[2])
  lims <- coord_cartesian(xlim = range(vec$PC1), ylim = range(vec$PC2))
  p <- ggplot(vec, aes(PC1, PC2, colour = factor(STATUS))) +
    geom_point(size = 2.6, alpha = .85) +
    scale_colour_manual(values = rep(PAL, 4), name = "Improvement status", na.value = "grey70") +
    labs(x = xl, y = yl) + lims + theme_pub(FS[["fig2"]])
  save_fig(p, "Fig2b_pca_PC1PC2_status", 7, 5.5)
  oil <- vec[!is.na(OIL_PCT)]
  p <- ggplot(oil, aes(PC1, PC2, colour = OIL_PCT)) +
    geom_point(size = 2.6) +
    scale_colour_viridis_c(option = "C", name = "Seed oil (%)") +
    labs(x = xl, y = yl) + lims + theme_pub(FS[["fig2"]])
  save_fig(p, "Fig2e_pca_by_oil_content", 7, 5.5)

  p <- ggplot(vec, aes(PC1, PC2, colour = factor(REGION))) +
    geom_point(size = 2.2, alpha = .85) +
    scale_colour_manual(values = rep(PAL, 4), name = "Region of origin", na.value = "grey70") +
    labs(x = xl, y = yl) + theme_pub(FS[["s9"]])
  save_fig(p, "Fig2b_pca_PC1PC2_region", 7, 5.5)

  ## Figure S1: 2-percent bins closed on the left, the mean as a dashed line and
  ## the summary statistics in a box
  o  <- oil$OIL_PCT
  st <- sprintf(paste("n = %d accessions", "Mean %.1f%% (SD %.1f)", "Median %.1f%%",
                      "Interquartile range %.0f to %.0f%%", "Range %.0f to %.0f%%", sep = "\n"),
                length(o), mean(o), sd(o), median(o), quantile(o, .25), quantile(o, .75), min(o), max(o))
  txt <- .9 * FS[["s1"]] / .pt
  s1 <- function(p) p +
    geom_vline(xintercept = mean(o), linewidth = .45, linetype = 2, colour = "grey15") +
    annotate("text", x = mean(o) + .4, y = Inf, label = sprintf("Mean %.1f%%", mean(o)),
             hjust = 0, vjust = 1.6, size = txt) +
    annotate("label", x = 42, y = Inf, label = st, hjust = 0, vjust = 1.08, size = txt,
             lineheight = 1.05, fill = "white") +
    scale_x_continuous(breaks = seq(15, 55, 5), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    coord_cartesian(xlim = c(15, 55)) +
    labs(x = "Seed oil content (%)", y = "Number of accessions") + theme_pub(FS[["s1"]])
  p <- s1(ggplot(oil, aes(OIL_PCT)) +
    geom_histogram(binwidth = 2, boundary = 0, closed = "left", fill = "#5B7DB1",
                   colour = "white", linewidth = .3))
  save_fig(p, "FigS8_oil_content_distribution", 6, 4)

  p <- ggplot(data.table(PC = seq_along(pct), V = pct)[1:15], aes(PC, V)) +
    geom_col(fill = PAL[1]) + geom_line(colour = PAL[4]) + geom_point(colour = PAL[4]) +
    scale_x_continuous(breaks = 1:15) +
    labs(x = "Principal component", y = "Variance explained (%)") + theme_pub(FS[["s6"]])
  save_fig(p, "Fig2a_pca_scree", 6, 4)
}

## ---- Figure S7: the BIONJ tree with the same key labels -----------------------
f_tr <- file.path(OUT, "06_tree", "sunflower_bionj.newick")
f_as <- file.path(AD, "group_assignment.txt")
if (have(f_tr) && have(f_as) && requireNamespace("ggtree", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(ape); library(ggtree) })
  tr  <- read.tree(f_tr)
  asg <- fread(f_as, header = FALSE, col.names = c("SAMPLE", "GROUP"))
  p <- ggtree(tr, layout = "fan", open.angle = 12, linewidth = .3)
  p <- p %<+% as.data.frame(asg) +
    geom_tippoint(aes(colour = GROUP), size = 1.6) +
    scale_cluster(name = "Cluster") +
    theme(legend.position = "right", legend.text = element_text(size = .8 * FS[["s7"]]),
          legend.title = element_text(size = FS[["s7"]]))
  save_fig(p, "Fig4_bionj_tree", 9, 9)
}

## ---- Figure 3: pairwise F_ST between clusters ---------------------------------
f_fst <- find_tab("pairwise_fst_clusters\\.csv$")
if (!file.exists(f_fst)) f_fst <- file.path(OUT, "04_diversity", "pairwise_fst_summary.tsv")
if (have(f_fst)) {
  fst <- fread(f_fst)
  fst[, WEIGHTED_FST := as.numeric(WEIGHTED_FST)]
  lv <- sort(unique(c(fst$POP1, fst$POP2)))
  fst[, `:=`(POP1 = factor(POP1, lv), POP2 = factor(POP2, lv))]
  p <- ggplot(fst, aes(POP1, POP2, fill = WEIGHTED_FST)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.3f", WEIGHTED_FST)), size = FS[["fig3"]] / .pt) +
    scale_fill_gradient(low = "#F7FBFF", high = PAL[4], name = expression(F[ST])) +
    labs(x = NULL, y = NULL) + theme_pub(FS[["fig3"]]) + theme(panel.grid.major = element_blank())
  save_fig(p, "Fig6_pairwise_fst", 7, 6)
}

## ---- Figure 4: F_ROH by cluster -----------------------------------------------
f_roh <- find_tab("_ROH\\.csv$")
if (have(f_roh)) {
  roh <- fread(f_roh)
  roh[, GROUP := factor(GROUP, levels = GRP_LEVELS)]
  p <- ggplot(roh, aes(GROUP, F_ROH, fill = GROUP)) +
    geom_boxplot(outlier.size = .8, alpha = .85) +
    scale_cluster("fill", guide = "none") +
    scale_x_discrete(labels = grp_label) +
    labs(x = NULL, y = expression(F[ROH])) + theme_pub(FS[["fig4"]])
  save_fig(p, "Fig7_FROH_by_group", 7, 4.5)
}

## ---- Figure 5: linkage disequilibrium decay -----------------------------------
f_ld <- find_tab("LD_binned_empirical\\.csv$")
if (have(f_ld)) {
  emp <- fread(f_ld)
  p <- ggplot(emp[DISTANCE_KB <= 500], aes(DISTANCE_KB, R2, colour = SET)) +
    geom_line(linewidth = .7) +
    geom_hline(yintercept = 0.2, linetype = 2, colour = "grey40") +
    scale_colour_manual(values = CLUSTER_COLS, name = NULL) +
    labs(x = "Physical distance (kb)",
         y = expression("Linkage disequilibrium," ~ r^2)) +
    theme_pub(FS[["fig5"]])
  save_fig(p, "Fig9_ld_decay", 7.5, 5)
  p2 <- p + scale_x_log10() + labs(x = "Physical distance (kb, log scale)")
  save_fig(p2, "Fig9b_ld_decay_log", 7.5, 5)
}

## ---- Figure S2: kinship ---------------------------------------------------------
f_grm <- file.path(QC, "grm.grm")
if (have(f_grm)) {
  gl <- fread(f_grm, col.names = c("i", "j", "nsnp", "k"))
  n  <- max(gl$i)
  K  <- matrix(0, n, n)
  K[cbind(gl$i, gl$j)] <- gl$k
  K[cbind(gl$j, gl$i)] <- gl$k
  gid <- fread(paste0(f_grm, ".id"), header = FALSE)
  dimnames(K) <- list(gid$V2, gid$V2)
  ord <- hclust(as.dist(1 - K / max(K)), method = "average")$order
  kl <- as.data.table(reshape2::melt(K[ord, ord], varnames = c("A", "B"), value.name = "K"))
  p <- ggplot(kl, aes(A, B, fill = pmin(K, quantile(K, .995)))) + geom_raster() +
    scale_fill_viridis_c(option = "B", name = "Kinship") +
    labs(x = NULL, y = NULL) + theme_pub(FS[["s2"]]) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank())
  save_fig(p, "Fig11_kinship_heatmap", 7.5, 6.5)
}

## ---- Figure S10: nucleotide diversity along the genome ------------------------
f_pi <- find_tab("windowed_pi\\.csv$")
if (have(f_pi)) {
  pw <- fread(f_pi)
  pw[, CHRL := factor(CHR_LAB[as.character(CHROM)], levels = CHR_LAB)]
  p <- ggplot(pw, aes(BIN_START / 1e6, PI, colour = GROUP)) +
    geom_line(linewidth = .3, alpha = .8) +
    facet_wrap(~CHRL, scales = "free_x", ncol = 3) +
    scale_cluster(name = "Cluster") +
    scale_y_continuous(labels = function(x) format(x * 1e5, trim = TRUE)) +
    labs(x = "Position (Mb)",
         y = expression(pi %*% 10^5 ~ "(1-Mb windows, 200-kb step)")) +
    theme_pub(FS[["s10"]]) +
    theme(legend.position = "bottom", panel.spacing.x = unit(14, "pt")) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))
  save_fig(p, "Fig8_pi_along_genome", 11, 12)
}

cat("figures written to", FIG, "\n")
