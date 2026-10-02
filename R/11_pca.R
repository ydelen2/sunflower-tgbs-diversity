## Step 11 - principal component analysis ---------------------------------------
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))

PC  <- file.path(OUT, "02_pca")
vec <- fread(file.path(PC, "pca.eigenvec"))
setnames(vec, 1:2, c("FID", "SAMPLE"))
val <- fread(file.path(PC, "pca.eigenval"), header = FALSE)$V1

## Variance explained is eigenvalue / trace(GRM), and the GRM must be built on
## the same markers as the PCA. plink2 --pca 20 writes only 20 eigenvalues, so
## sum(val) is not the total.
grm_f <- file.path(PC, "pca_grm.grm")
if (!file.exists(grm_f)) grm_f <- file.path(OUT, "01_qc", "grm.grm")
if (file.exists(grm_f)) {
  gl    <- fread(grm_f, col.names = c("i", "j", "nsnp", "k"))
  trace <- gl[i == j, sum(k)]
  pct   <- round(100 * val / trace, 2)
  cat(sprintf("GRM trace = %.2f over %d accessions; PC1-3 explain %.2f/%.2f/%.2f%% of total variance\n",
              trace, gl[i == j, .N], pct[1], pct[2], pct[3]))
  cat(sprintf("(the top-20 PCs together account for %.1f%% of total variance)\n",
              sum(pct[seq_len(min(20, length(pct)))])))
} else {
  warning("grm.grm not found - falling back to the top-20 denominator; ",
          "axis labels will overstate variance explained")
  pct <- round(100 * val / sum(val), 2)
}

meta <- read_meta()
if (!is.null(meta)) vec <- merge(vec, meta, by = "SAMPLE", all.x = TRUE)

## scree
p0 <- ggplot(data.table(PC = seq_along(pct), V = pct)[1:15], aes(PC, V)) +
  geom_col(fill = PAL[1]) + geom_line(colour = PAL[4]) + geom_point(colour = PAL[4]) +
  scale_x_continuous(breaks = 1:15) +
  labs(x = "Principal component", y = "Variance explained (%)") + theme_pub()
save_fig(p0, "Fig2a_pca_scree", 6, 4)

## legend titles for the reader, not the column names
GRP_LAB <- c(STATUS = "Improvement status", REGION = "Region of origin",
             COUNTRY = "Country of origin")

plot_pc <- function(dt, x, y, colvar) {
  lab <- if (is.null(colvar)) NULL else
    if (!is.na(GRP_LAB[colvar])) unname(GRP_LAB[colvar]) else colvar
  ggplot(dt, aes(.data[[x]], .data[[y]],
                 colour = if (is.null(colvar)) NULL else factor(.data[[colvar]]))) +
    geom_point(size = 2.2, alpha = .85) +
    scale_colour_manual(values = rep(PAL, 4), name = lab, na.value = "grey70") +
    labs(x = sprintf("%s (%.2f%%)", x, pct[as.integer(sub("PC", "", x))]),
         y = sprintf("%s (%.2f%%)", y, pct[as.integer(sub("PC", "", y))])) +
    theme_pub()
}

## color by whatever grouping columns exist in the metadata
grpcols <- intersect(c("STATUS", "REGION", "COUNTRY"), names(vec))
if (length(grpcols) == 0) grpcols <- NA_character_

for (g in grpcols) {
  sfx <- if (is.na(g)) "nogroup" else tolower(g)
  save_fig(plot_pc(vec, "PC1", "PC2", if (is.na(g)) NULL else g), paste0("Fig2b_pca_PC1PC2_", sfx), 7, 5.5)
  save_fig(plot_pc(vec, "PC1", "PC3", if (is.na(g)) NULL else g), paste0("Fig2c_pca_PC1PC3_", sfx), 7, 5.5)
  save_fig(plot_pc(vec, "PC2", "PC3", if (is.na(g)) NULL else g), paste0("Fig2d_pca_PC2PC3_", sfx), 7, 5.5)
}

## Seed oil content as a continuous covariate.
if (!is.null(meta) && "OIL_PCT" %in% names(vec)) {
  vec[, OIL_PCT := suppressWarnings(as.numeric(OIL_PCT))]
  oil <- vec[!is.na(OIL_PCT)]
  cat(sprintf("\nOil content available for %d of %d genotyped accessions\n",
              nrow(oil), nrow(vec)))
  for (k in 1:5) {
    ct <- cor.test(oil[[paste0("PC", k)]], oil$OIL_PCT)
    cat(sprintf("  PC%d ~ oil%%: r = %+.3f  p = %.3g\n", k, ct$estimate, ct$p.value))
  }
  pO <- ggplot(oil, aes(PC1, PC2, colour = OIL_PCT)) +
    geom_point(size = 2.4) +
    scale_colour_viridis_c(option = "C", name = "Seed oil (%)") +
    labs(x = sprintf("PC1 (%.2f%%)", pct[1]), y = sprintf("PC2 (%.2f%%)", pct[2])) +
    theme_pub()
  save_fig(pO, "Fig2e_pca_by_oil_content", 7, 5.5)

  pH <- ggplot(oil, aes(OIL_PCT)) +
    geom_histogram(binwidth = 2, fill = PAL[3], colour = "white") +
    labs(x = "Seed oil content (%)", y = "Number of accessions",
         subtitle = NULL) +
    theme_pub()
  save_fig(pH, "FigS8_oil_content_distribution", 6, 4)
}

## STATUS x REGION cross-tabulation (non-orthogonal in this panel).
if (!is.null(meta) && all(c("STATUS", "REGION") %in% names(vec))) {
  cat("\nSTATUS x REGION cross-tabulation:\n")
  print(table(vec$STATUS, vec$REGION))
}

fwrite(vec[, !"FID"], file.path(TAB, "Data4_pca_scores.csv"))
fwrite(data.table(PC = seq_along(pct), EIGENVALUE = val, VAR_PCT = pct),
       file.path(TAB, "TableS8_pca_variance.csv"))

## ANOVA of PC1 on each a priori grouping
if (!is.null(meta)) {
  for (g in grpcols[!is.na(grpcols)]) {
    f <- as.formula(sprintf("PC1 ~ factor(%s)", g))
    a <- summary(aov(f, data = vec[!is.na(get(g))]))
    cat("\n== PC1 ~", g, "==\n"); print(a)
  }
}
cat("step 11 done\n")
