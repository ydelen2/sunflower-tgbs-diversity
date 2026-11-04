## Step 26 - Figure S10: linkage disequilibrium decay for each chromosome -------
## (a) the whole-panel decay curve of every chromosome, from the per-chromosome
##     runs of slurm/12_ld_by_chromosome.slurm, drawn as the isotonic fit to the
##     2-kb bin means from which Table S14 is read, with the whole-panel fit of
##     Figure 4 for reference; colours rank the chromosomes by r-squared at 100 kb
## (b) mean r-squared at 100 kb per chromosome, read from an isotonic fit to the
##     2-kb bin means as in Table S14, against the whole-panel value
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
library(grid)

LDC <- file.path(OUT, "05_ld", "bychr")
files <- list.files(LDC, pattern = "^chr_Ha[0-9]+\\.stat(\\.gz)?$", full.names = TRUE)
stopifnot(length(files) == 17L)

read_ld <- function(f) {
  d <- fread(cmd = if (grepl("\\.gz$", f)) paste("zcat", shQuote(f)) else NULL,
             file = if (grepl("\\.gz$", f)) NULL else f)
  setnames(d, 1:2, c("Dist", "MeanR2"))
  np <- if (ncol(d) >= 6) d[[6]] else rep(1L, nrow(d))
  data.table(dist_bp = as.numeric(d$Dist), r2 = as.numeric(d$MeanR2), n = as.numeric(np))
}

BIN <- 2000
ld <- rbindlist(lapply(files, function(f) {
  d <- read_ld(f)
  d[, CHROMOSOME := sub("\\.stat(\\.gz)?$", "", sub("^chr_", "", basename(f)))]
  d
}))
ld <- ld[dist_bp > 0 & is.finite(r2)]
emp <- ld[, .(R2 = sum(r2 * n) / sum(n), PAIRS = sum(n)),
          by = .(CHROMOSOME, bin = floor(dist_bp / BIN))]
emp[, DISTANCE_KB := (bin + 0.5) * BIN / 1000]
setorder(emp, CHROMOSOME, bin)

## the monotone (isotonic) fit of every chromosome, and r-squared at 100 kb from it, as in Table S14
fit <- emp[, {
  o <- order(bin)
  .(DISTANCE_KB = DISTANCE_KB[o], R2 = -isoreg(bin[o], -R2[o])$yf)
}, by = CHROMOSOME]
at100 <- fit[, .(R2_100KB = R2[which.min(abs(DISTANCE_KB - 100))]), by = CHROMOSOME]
at100[, ORD := as.integer(sub("Ha", "", CHROMOSOME))]
setorder(at100, ORD)
chr_levels <- at100$CHROMOSOME
fit[, CHROMOSOME := factor(CHROMOSOME, levels = chr_levels)]
at100[, CHROMOSOME := factor(CHROMOSOME, levels = chr_levels)]

## the whole-panel curve and its value at 100 kb (Data file 9, Table 2)
f_ld <- list.files(TAB, "LD_binned_empirical\\.csv$", full.names = TRUE)[1]
all_curve <- NULL; all100 <- NA_real_
if (!is.na(f_ld) && file.exists(f_ld)) {
  e <- fread(f_ld)[SET == "All accessions"]
  o <- order(e$DISTANCE_KB); mono <- -isoreg(e$DISTANCE_KB[o], -e$R2[o])$yf
  all_curve <- data.table(DISTANCE_KB = e$DISTANCE_KB[o], R2 = mono)
  all100 <- mono[which.min(abs(e$DISTANCE_KB[o] - 100))]
}

## one colour per chromosome, ordered by r-squared at 100 kb (highest first); the legend lists that order
rank_levels <- as.character(at100[order(-R2_100KB), CHROMOSOME])
chr_cols <- setNames(scales::viridis_pal(option = "H")(17), rank_levels)
pa <- ggplot(fit[DISTANCE_KB <= 500], aes(DISTANCE_KB, R2, colour = CHROMOSOME)) +
  geom_line(linewidth = .7) +
  scale_colour_manual(values = chr_cols, breaks = rank_levels, name = NULL) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(linewidth = 1.6))) +
  labs(x = "Physical distance (kb)", y = expression("Linkage disequilibrium," ~ r^2)) +
  theme_pub(12) + theme(legend.key.height = unit(.42, "cm"))
if (!is.null(all_curve)) {
  pa <- pa + geom_line(data = all_curve[DISTANCE_KB <= 500], aes(DISTANCE_KB, R2),
                       inherit.aes = FALSE, colour = "black", linewidth = .9, linetype = 1)
}

pb <- ggplot(at100, aes(CHROMOSOME, R2_100KB)) +
  geom_segment(aes(xend = CHROMOSOME, y = 0, yend = R2_100KB), colour = "grey70", linewidth = .6) +
  geom_point(aes(colour = CHROMOSOME), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = chr_cols) +
  labs(x = NULL, y = expression("Mean" ~ r^2 ~ "at 100 kb")) +
  theme_pub(12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (is.finite(all100)) pb <- pb + geom_hline(yintercept = all100, linetype = 2, colour = "grey40")

draw <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1.15, 1), "null"))))
  print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pb, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.text("(a)", x = unit(0.012, "npc"), y = unit(0.985, "npc"), just = c("left", "top"),
            gp = gpar(fontface = "bold", fontsize = 13))
  grid.text("(b)", x = unit(0.012, "npc"), y = unit(0.455, "npc"), just = c("left", "top"),
            gp = gpar(fontface = "bold", fontsize = 13))
}
cairo_pdf(file.path(FIG, "FigS10_ld_by_chromosome.pdf"), width = 7.5, height = 8.5); draw(); dev.off()
png(file.path(FIG, "FigS10_ld_by_chromosome.png"), width = 7.5, height = 8.5, units = "in", res = 400, type = "cairo"); draw(); dev.off()
print(at100[, .(CHROMOSOME, R2_100KB = round(R2_100KB, 3))])
cat("whole panel at 100 kb:", round(all100, 3), "\n")
cat("Figure S10 written to", FIG, "\n")
