## Masking experiment on imputation (slurm/15): accuracy of the re-imputed
## genotypes at the masked cells, and genomic inbreeding from runs of
## homozygosity on each re-imputed set against the provider's imputation.
suppressPackageStartupMessages(library(data.table))
PROJ <- Sys.getenv("PROJ", "/work/dweikat/ydelen2/sunflower_tgbs")
ANA <- file.path(PROJ, "analysis"); OUT <- file.path(PROJ, "results")
MK <- file.path(OUT, "12_review", "masking"); TAB <- file.path(OUT, "09_tables")
GENOME_KB <- sum(fread(file.path(ANA, "ref/chr_lengths.txt"), header = FALSE)$V2) / 1000
CL <- c("Q1", "Q2", "Q3", "Q4", "Q5")

grp <- fread(file.path(OUT, "03_admixture", "group_assignment.txt"), header = FALSE, col.names = c("SAMPLE", "GROUP"))
samples <- readLines(file.path(MK, "samples.txt"))
froh <- function(f) { h <- fread(f); h[, .(SAMPLE = IID, F_ROH = KB / GENOME_KB, NSEG, KBAVG)] }
ref <- froh(file.path(OUT, "04_diversity", "roh.hom.indiv"))

gz <- function(f) fread(cmd = paste("gzip -dc", shQuote(f)))              # no R.utils needed
code_gt <- function(g) { g <- gsub("|", "/", g, fixed = TRUE); c("0/0" = 0L, "0/1" = 1L, "1/0" = 1L, "1/1" = 2L)[g] }

levels_ <- c(mask0 = 0, mask10 = 0.10, mask20 = 0.20)
summ <- list(); bycl <- list(); byacc <- list()
for (tag in names(levels_)) {
  fr <- froh(file.path(MK, sprintf("roh_%s.hom.indiv", tag)))
  m <- merge(merge(fr, ref[, .(SAMPLE, F_ROH_PROVIDER = F_ROH)], by = "SAMPLE"), grp, by = "SAMPLE")
  m[, LEVEL := tag]
  byacc[[tag]] <- m
  cm <- m[, .(N = .N, F_ROH = mean(F_ROH), F_ROH_PROVIDER = mean(F_ROH_PROVIDER)), by = GROUP]
  acc_all <- acc_hom <- acc_het <- NA_real_; n_masked <- 0L; cl_acc <- NULL
  if (levels_[[tag]] > 0) {
    masked <- gz(file.path(MK, sprintf("%s.masked.tsv.gz", tag)))
    gt <- gz(file.path(MK, sprintf("%s.imputed.gt.tsv.gz", tag)))          # CHROM, POS, one integer dosage per sample
    setnames(gt, c("CHROM", "POS", samples))
    key <- paste(gt$CHROM, gt$POS)
    masked[, ROW := match(paste(CHROM, POS), key)]
    stopifnot(!anyNA(masked$ROW))
    imp <- integer(nrow(masked))
    for (s in unique(masked$SAMPLE)) { i <- which(masked$SAMPLE == s); imp[i] <- gt[[s]][masked$ROW[i]] }
    masked[, IMPUTED := imp]
    masked[, GT := code_gt(GT)]
    masked[, OK := GT == IMPUTED]
    masked[, HET := GT == 1L]
    masked <- merge(masked, grp, by = "SAMPLE")
    n_masked <- nrow(masked)
    acc_all <- masked[, mean(OK)]; acc_hom <- masked[HET == FALSE, mean(OK)]; acc_het <- masked[HET == TRUE, mean(OK)]
    cl_acc <- masked[, .(ACCURACY = mean(OK), ACCURACY_HOM = mean(OK[!HET]), ACCURACY_HET = mean(OK[HET]),
                         SHARE_HET_MASKED = mean(HET)), by = GROUP]
    rm(gt, masked); gc()
  }
  summ[[tag]] <- data.table(LEVEL = tag, FRACTION_MASKED = levels_[[tag]], N_MASKED = n_masked,
                            ACCURACY = acc_all, ACCURACY_HOMOZYGOUS = acc_hom, ACCURACY_HETEROZYGOUS = acc_het,
                            MEAN_F_ROH = mean(m$F_ROH), MEAN_F_ROH_PROVIDER = mean(m$F_ROH_PROVIDER),
                            PEARSON_WITH_PROVIDER = cor(m$F_ROH, m$F_ROH_PROVIDER),
                            SPEARMAN_WITH_PROVIDER = cor(m$F_ROH, m$F_ROH_PROVIDER, method = "spearman"),
                            CLUSTER_ORDER = paste(cm[GROUP != "admixed"][order(-F_ROH), GROUP], collapse = " > "),
                            CLUSTER_ORDER_PROVIDER = paste(cm[GROUP != "admixed"][order(-F_ROH_PROVIDER), GROUP], collapse = " > "))
  cc <- if (is.null(cl_acc)) cm else merge(cm, cl_acc, by = "GROUP")
  cc[, LEVEL := tag]; bycl[[tag]] <- cc
}
summ <- rbindlist(summ); bycl <- rbindlist(bycl, fill = TRUE); byacc <- rbindlist(byacc)
setcolorder(bycl, c("LEVEL", "GROUP"))
print(summ); print(bycl)
fwrite(summ[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], file.path(TAB, "review_masking.csv"))
fwrite(bycl[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], file.path(TAB, "review_masking_by_cluster.csv"))
fwrite(byacc[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], file.path(TAB, "review_masking_by_accession.csv"))
cat("step 25 done\n")
