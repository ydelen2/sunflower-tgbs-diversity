## Step 10 - accession- and site-level QC ---------------------------------------
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))

QC <- file.path(OUT, "01_qc")

## sample_stats.tsv        - step 01, all 276 delivered columns, no keep list
## sample_stats_kept.tsv   - step 02, the analysis accessions only
ss_all <- fread(file.path(QC, "sample_stats.tsv"))
kf <- file.path(QC, "sample_stats_kept.tsv")
ss <- if (file.exists(kf)) fread(kf) else ss_all
cat(sprintf("accessions: %d delivered, %d in the analysis panel\n", nrow(ss_all), nrow(ss)))

## Ho is heterozygous/CALLED. heterozygous/ALL-SITES equals Ho * (1 - missing)
## and therefore correlates with missingness by construction; both are computed
## so the Methods can state the difference.
ss[, HET_PER_ALL := HO * (1 - MISS_RATE)]
for (v in c("HO", "HET_PER_ALL", "CORE_HO")) {
  ct <- cor.test(ss$MISS_RATE, ss[[v]])
  cat(sprintf("%-12s ~ missingness : r = %+.3f (p = %.3g)\n", v, ct$estimate, ct$p.value))
}

pl <- melt(ss[, .(SAMPLE, MISS_RATE,
                  `Heterozygous / called sites` = HO,
                  `Heterozygous / called, common sites` = CORE_HO,
                  `Heterozygous / all sites` = HET_PER_ALL)],
           id.vars = c("SAMPLE", "MISS_RATE"), variable.name = "which", value.name = "Ho")
p1 <- ggplot(pl, aes(MISS_RATE, Ho)) +
  geom_point(alpha = .6, size = 1.6, colour = PAL[1]) +
  geom_smooth(method = "lm", se = TRUE, colour = PAL[4], linewidth = .6) +
  facet_wrap(~which) +
  labs(x = "Missing genotype rate per accession", y = "Observed heterozygosity") +
  theme_pub()
save_fig(p1, "FigS1_het_vs_missing", 10, 4)

thr_samp <- suppressWarnings(as.numeric(Sys.getenv("MAX_SAMPLE_MISS", "")))
p2 <- ggplot(ss_all, aes(MISS_RATE)) +
  geom_histogram(bins = 40, fill = PAL[1], colour = "white") +
  labs(x = "Missing genotype rate per accession", y = "Number of accessions") + theme_pub()
if (is.finite(thr_samp))
  p2 <- p2 + geom_vline(xintercept = thr_samp, linetype = 2, colour = PAL[4])
save_fig(p2, "FigS2_sample_missingness", 6, 4)

## --- site level ---------------------------------------------------------------
## site_stats_filtered.tsv is written by step 02 on the analysis accessions and
## carries a PASS flag; every site-level statistic in the manuscript describes
## the PASS subset, which is the marker set the analyses actually used.
sf <- file.path(QC, "site_stats_filtered.tsv")
stopifnot(file.exists(sf))
st_all <- fread(sf)
stopifnot("PASS" %in% names(st_all))
st <- st_all[PASS == 1]
st[, CHRL := factor(CHR_LAB[as.character(CHROM)], levels = CHR_LAB)]
cat(sprintf("sites: %d candidates, %d retained\n", nrow(st_all), nrow(st)))

## sites removed by each criterion, and the retained set, with their means
pic <- function(p) 1 - (p^2 + (1 - p)^2) - 2 * p^2 * (1 - p)^2
st_all[, `:=`(CLASS = fifelse(PASS == 1, "retained", as.character(REASON)), PIC = pic(MAF))]
st[, PIC := pic(MAF)]
cls <- st_all[, .(N_SITES = .N, MAF_mean = round(mean(MAF), 4), He_mean = round(mean(HE), 4),
                  Ho_mean = round(mean(HO), 4), PIC_mean = round(mean(PIC), 4)),
              by = CLASS][order(-N_SITES)]
print(cls)
fwrite(cls, file.path(TAB, "TableS3_site_filter_classes.csv"))
cr <- file.path(QC, "control_replicate_concordance.tsv")
if (file.exists(cr)) fwrite(fread(cr), file.path(TAB, "replicate_concordance.csv"))

p3 <- ggplot(st, aes(MAF)) + geom_histogram(binwidth = .01, fill = PAL[3], colour = "white") +
  labs(x = "Minor allele frequency", y = "Number of SNPs") + theme_pub()
p4 <- ggplot(st, aes(MISS)) + geom_histogram(binwidth = .02, fill = PAL[2], colour = "white") +
  labs(x = "Missing rate per SNP", y = "Number of SNPs") + theme_pub()
p5 <- ggplot(st, aes(HO)) + geom_histogram(binwidth = .01, fill = PAL[5], colour = "white") +
  labs(x = "Heterozygous call rate per SNP", y = "Number of SNPs") + theme_pub()
save_fig(p3, "FigS3_maf_spectrum"); save_fig(p4, "FigS4_site_missing"); save_fig(p5, "FigS5_site_het")

dens <- st[!is.na(CHRL), .N, by = .(CHRL, WIN = floor(POS / 1e6))]
p6 <- ggplot(dens, aes(WIN, CHRL, fill = N)) + geom_tile(height = .8) +
  scale_fill_viridis_c(name = "SNP / Mb", option = "C") +
  scale_y_discrete(limits = rev(levels(dens$CHRL))) +
  labs(x = "Position (Mb)", y = NULL) + theme_pub() +
  theme(panel.grid = element_blank())
save_fig(p6, "Fig1_snp_density", 9, 5)

## --- Table 1: marker summary, retained markers only ---------------------------
len <- fread(file.path(ANA, "ref/chr_lengths.txt"), col.names = c("CHROM", "LEN"))
len[, CHRL := factor(CHR_LAB[as.character(CHROM)], levels = CHR_LAB)]

tab1 <- st[, .(N_SNP = .N, MAF_mean = mean(MAF), MAF_median = median(MAF),
               He_mean = mean(HE), Ho_mean = mean(HO), PIC_mean = mean(PIC),
               Missing_mean = mean(MISS)), by = CHRL][order(CHRL)]
tab1 <- merge(tab1, len[, .(CHRL, LEN_Mb = LEN / 1e6)], by = "CHRL", all.x = TRUE)
tab1[, `:=`(SNP_per_Mb = N_SNP / LEN_Mb, mean_spacing_kb = LEN_Mb * 1000 / N_SNP)]
setcolorder(tab1, c("CHRL", "LEN_Mb", "N_SNP", "SNP_per_Mb", "mean_spacing_kb"))

total <- st[, .(CHRL = "Total", LEN_Mb = sum(len$LEN) / 1e6, N_SNP = .N,
                SNP_per_Mb = .N / (sum(len$LEN) / 1e6),
                mean_spacing_kb = (sum(len$LEN) / 1000) / .N,
                MAF_mean = mean(MAF), MAF_median = median(MAF),
                He_mean = mean(HE), Ho_mean = mean(HO), PIC_mean = mean(PIC),
                Missing_mean = mean(MISS))]
tab1 <- rbind(tab1, total, fill = TRUE)
fwrite(tab1, file.path(TAB, "Table1_marker_summary.csv"))
print(tab1)

fwrite(ss, file.path(TAB, "Data2_accession_qc.csv"))
cat("step 10 done\n")
