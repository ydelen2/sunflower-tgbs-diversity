## Step 14 - diversity indices, AMOVA, Fst tables, ROH, pi ----------------------
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
suppressPackageStartupMessages({
  library(hierfstat); library(poppr); library(adegenet); library(vcfR)
})

DV   <- file.path(OUT, "04_diversity")
asg  <- fread(file.path(OUT, "03_admixture/group_assignment.txt"),
              header = FALSE, col.names = c("SAMPLE", "GROUP"))
meta <- read_meta()
if (!is.null(meta)) asg <- merge(asg, meta, by = "SAMPLE", all.x = TRUE)
grps <- sort(unique(asg$GROUP))

## ---------------------------------------------------------------------------
## A. Per-group diversity from the vcftools output (all filtered SNPs)
## ---------------------------------------------------------------------------
read_frq <- function(f) {
  fr <- fread(f, skip = 1, col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "P1", "P2"))
  ## A site with no called genotype in the group has no defined frequency.
  ## Keeping such rows makes every frequency-derived mean NA in the small
  ## clusters and understates the polymorphic fraction, so they are dropped and
  ## the statistics are computed over the sites the group actually carries.
  fr <- fr[N_CHR > 0 & is.finite(P1) & is.finite(P2)]
  fr[, `:=`(He  = 2 * P1 * P2,
            MAF = pmin(P1, P2),
            PIC = 1 - (P1^2 + P2^2) - 2 * P1^2 * P2^2)]   # Botstein PIC, biallelic
  fr[]
}

div <- rbindlist(lapply(grps, function(g) {
  ff <- file.path(DV, sprintf("grp_%s.frq", g))
  fh <- file.path(DV, sprintf("grp_%s.het", g))
  if (!file.exists(ff) || !file.exists(fh)) return(NULL)
  fr <- read_frq(ff)
  h  <- fread(fh)
  setnames(h, gsub("[()]", "", names(h)))
  data.table(GROUP         = g,
             N             = nrow(h),
             N_SNP_called  = nrow(fr),
             N_polymorphic = fr[MAF > 0, .N],
             PctPoly       = 100 * fr[MAF > 0, .N] / nrow(fr),
             MAF           = mean(fr$MAF),
             He            = mean(fr$He),
             PIC           = mean(fr$PIC),
             Ho            = mean(1 - h$OHOM / h$N_SITES),
             Fis_vcftools  = mean(h$F))
}))
## No 1 - Ho/He column here. Expected heterozygosity is averaged over every
## site called in the cluster, including those fixed within it, while observed
## heterozygosity comes from the per-individual vcftools calculation over each
## individual's own called sites. The two denominators differ, so their ratio is
## not a valid inbreeding coefficient; the vcftools F, which is internally
## consistent, is reported instead.

## Cluster identity: seed oil content and the composition of each cluster in
## terms of the a priori categories. This is what makes the clusters
## interpretable, so it belongs in the same table as their diversity.
if (!is.null(meta)) {
  m2 <- merge(asg[, .(SAMPLE, GROUP)], meta, by = "SAMPLE", all.x = TRUE)
  m2[, OIL := suppressWarnings(as.numeric(OIL_PCT))]
  ident <- m2[, .(
    N_oil        = sum(!is.na(OIL)),
    oil_mean     = round(mean(OIL, na.rm = TRUE), 1),
    oil_sd       = round(sd(OIL, na.rm = TRUE), 1),
    n_oilseed    = sum(TYPE == "oilseed", na.rm = TRUE),
    n_confection = sum(TYPE == "confection", na.rm = TRUE),
    n_restorer   = sum(HETEROTIC_GROUP == "R", na.rm = TRUE),
    n_maintainer = sum(HETEROTIC_GROUP == "B", na.rm = TRUE),
    n_breeding   = sum(STATUS == "Breeding line", na.rm = TRUE),
    n_cultivar   = sum(STATUS == "Cultivar", na.rm = TRUE),
    n_landrace   = sum(STATUS == "Landrace", na.rm = TRUE),
    top_region   = names(sort(table(REGION), decreasing = TRUE))[1],
    top_region_n = max(table(REGION))
  ), by = GROUP]
  ## the admixed accessions have no per-group diversity files; they enter the
  ## table for their composition only
  div <- merge(div, ident, by = "GROUP", all = TRUE)
  div[is.na(N), N := asg[GROUP == .BY$GROUP, .N], by = GROUP]

  cat("\nseed oil content by cluster\n")
  print(summary(aov(OIL ~ factor(GROUP), data = m2[!is.na(OIL)])))
}

fwrite(div, file.path(TAB, "Table1_cluster_diversity_and_identity.csv"))
print(div)

## whole-collection row
fr_all <- read_frq(file.path(DV, "all.frq"))
cat(sprintf("\nWhole collection: %d SNPs | MAF %.3f | He %.3f | PIC %.3f\n",
            nrow(fr_all), mean(fr_all$MAF), mean(fr_all$He), mean(fr_all$PIC)))

## ---------------------------------------------------------------------------
## B. Pairwise Fst heat map
## ---------------------------------------------------------------------------
fstf <- file.path(DV, "pairwise_fst_summary.tsv")
if (file.exists(fstf)) {
  fst <- fread(fstf)
  fst[, WEIGHTED_FST := as.numeric(WEIGHTED_FST)]
  lv <- sort(unique(c(fst$POP1, fst$POP2)))
  fst[, `:=`(POP1 = factor(POP1, lv), POP2 = factor(POP2, lv))]
  p <- ggplot(fst, aes(POP1, POP2, fill = WEIGHTED_FST)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.3f", WEIGHTED_FST)), size = 3) +
    scale_fill_gradient(low = "#F7FBFF", high = PAL[4], name = expression(F[ST])) +
    labs(x = NULL, y = NULL) + theme_pub() + theme(panel.grid = element_blank())
  save_fig(p, "Fig6_pairwise_fst", 7, 6)
  fwrite(fst, file.path(TAB, "pairwise_fst_clusters.csv"))
}

## ---------------------------------------------------------------------------
## C. AMOVA on a random SNP subset (poppr cannot take the full marker set)
## ---------------------------------------------------------------------------
vcf  <- read.vcfR(file.path(OUT, "00_vcf/sunflower.filtered.vcf.gz"), verbose = FALSE)
set.seed(42)
keep <- sort(sample(nrow(vcf@fix), min(20000L, nrow(vcf@fix))))
gid  <- vcfR2genind(vcf[keep, ])
rm(vcf); gc()
## AMOVA is run on the assigned accessions only, matching the pairwise F_ST
## analysis. "admixed" is a residual category below the Q >= 0.60 threshold,
## not a deme, and including it as a group would test a different hypothesis.
gid_all <- gid            # every genotyped accession, used for the a priori
                          # groupings that do not depend on the clustering
pops <- asg[GROUP != "admixed"]
gid  <- gid[indNames(gid) %in% pops$SAMPLE]
cat(sprintf("AMOVA on %d assigned accessions in %d clusters\n",
            nrow(pops), uniqueN(pops$GROUP)))

strat_cols <- intersect(c("GROUP", "REGION", "COUNTRY", "STATUS"), names(pops))
st <- as.data.frame(pops[match(indNames(gid), SAMPLE), ..strat_cols])
st[] <- lapply(st, function(x) factor(ifelse(is.na(x), "unknown", as.character(x))))
strata(gid) <- st

## AMOVA results as a table rather than a printed R object
amova_tidy <- function(am, sig, model) {
  res <- as.data.frame(am$results)
  cvc <- as.data.frame(am$componentsofcovariance)
  phi <- as.data.frame(am$statphi)
  na <- function(n) rep(NA_real_, n)
  out <- rbind(
    data.table(MODEL = model, SECTION = "Variance", TERM = trimws(rownames(res)),
               DF = as.numeric(res[[1]]), SUM_SQ = round(as.numeric(res[[2]]), 3),
               MEAN_SQ = round(as.numeric(res[[3]]), 3), SIGMA = na(nrow(res)),
               PERCENT = na(nrow(res)), STATISTIC = na(nrow(res)),
               P_VALUE = na(nrow(res))),
    data.table(MODEL = model, SECTION = "Components of covariance",
               TERM = trimws(rownames(cvc)), DF = na(nrow(cvc)),
               SUM_SQ = na(nrow(cvc)), MEAN_SQ = na(nrow(cvc)),
               SIGMA = round(as.numeric(cvc[[1]]), 5),
               PERCENT = round(as.numeric(cvc[[2]]), 3),
               STATISTIC = na(nrow(cvc)), P_VALUE = na(nrow(cvc))),
    data.table(MODEL = model, SECTION = "Phi statistics",
               TERM = trimws(rownames(phi)), DF = na(nrow(phi)),
               SUM_SQ = na(nrow(phi)), MEAN_SQ = na(nrow(phi)),
               SIGMA = na(nrow(phi)), PERCENT = na(nrow(phi)),
               STATISTIC = round(as.numeric(phi[[1]]), 5), P_VALUE = na(nrow(phi))))
  if (!is.null(sig)) {
    pv <- as.numeric(sig$pvalue)
    nm <- if (!is.null(sig$names)) as.character(sig$names) else "Among groups"
    out <- rbind(out, data.table(MODEL = model, SECTION = "Permutation test",
      TERM = rep(nm, length.out = length(pv)), DF = na(length(pv)),
      SUM_SQ = na(length(pv)), MEAN_SQ = na(length(pv)), SIGMA = na(length(pv)),
      PERCENT = na(length(pv)),
      STATISTIC = round(as.numeric(sig$obs), 5)[seq_along(pv)], P_VALUE = pv))
  }
  out[]
}

gc_obj <- as.genclone(gid)
amres  <- poppr.amova(gc_obj, ~GROUP, within = FALSE, quiet = TRUE)
print(amres)
set.seed(1)
amsig <- randtest(amres, nrepet = 999)
print(amsig)

## hierarchical AMOVA when passport data give a second level
if ("REGION" %in% strat_cols && nlevels(st$REGION) > 1) {
  am2 <- poppr.amova(gc_obj, ~REGION/GROUP, within = FALSE, quiet = TRUE)
  print(am2)
} else am2 <- NULL
fwrite(rbind(amova_tidy(amres, amsig, "Ancestry clusters"),
             if (!is.null(am2)) amova_tidy(am2, NULL, "Region / clusters")),
       file.path(TAB, "TableS10c_AMOVA.csv"))

## hierfstat basic stats on the same subset
hf <- genind2hierfstat(gid, pop = factor(pops[match(indNames(gid), SAMPLE), GROUP]))
bstat <- basic.stats(hf)
print(bstat$overall)
fwrite(as.data.table(as.list(bstat$overall)),
       file.path(TAB, "TableS10_hierfstat_overall.csv"))

## Pairwise F_ST is taken from vcftools on the full marker set (Table S10).
## hierfstat::pairwise.WCfst on the AMOVA subset would repeat the same
## statistic with fewer markers at a large time cost, so it is not run here.

## ---------------------------------------------------------------------------
## D. Runs of homozygosity -> F_ROH
## ---------------------------------------------------------------------------
rf <- file.path(DV, "roh.hom.indiv")
if (file.exists(rf)) {
  roh <- fread(rf)
  GENOME_KB <- sum(fread(file.path(ANA, "ref/chr_lengths.txt"))$V2) / 1000
  roh[, F_ROH := KB / GENOME_KB]
  roh <- merge(roh[, .(SAMPLE = IID, NSEG, KB, KBAVG, F_ROH)], asg, by = "SAMPLE")
  p2 <- ggplot(roh, aes(GROUP, F_ROH, fill = GROUP)) +
    geom_boxplot(outlier.size = .8, alpha = .85) +
    scale_fill_manual(values = CLUSTER_COLS, guide = "none") +
    labs(x = NULL, y = expression(F[ROH])) + theme_pub()
  save_fig(p2, "Fig7_FROH_by_group", 7, 4.5)
  fwrite(roh, file.path(TAB, "Data7_ROH.csv"))
  print(roh[, .(n = .N, meanF = mean(F_ROH), meanNSEG = mean(NSEG)), by = GROUP][order(-meanF)])
}

## ---------------------------------------------------------------------------
## E. Nucleotide diversity along the genome
##    vcftools --window-pi divides by full window length; tGBS interrogates
##    ~1 site per 12 kb. Within-dataset comparison only.
## ---------------------------------------------------------------------------
pw <- rbindlist(lapply(grps, function(g) {
  f <- file.path(DV, sprintf("grp_%s.windowed.pi", g))
  if (!file.exists(f)) return(NULL)
  d <- fread(f); d[, GROUP := g]; d
}))
if (nrow(pw)) {
  pw[, CHRL := factor(CHR_LAB[as.character(CHROM)], levels = CHR_LAB)]
  p3 <- ggplot(pw, aes(BIN_START / 1e6, PI, colour = GROUP)) +
    geom_line(linewidth = .3, alpha = .8) +
    facet_wrap(~CHRL, scales = "free_x", ncol = 3) +
    scale_colour_manual(values = CLUSTER_COLS, name = "Cluster") +
    scale_y_continuous(labels = function(x) format(x * 1e5, trim = TRUE)) +
    labs(x = "Position (Mb)",
         y = expression(pi %*% 10^5 ~ "(1-Mb windows, 200-kb step)")) + theme_pub(9)
  save_fig(p3, "Fig8_pi_along_genome", 11, 11)
  print(pw[, .(mean_pi = mean(PI), median_pi = median(PI)), by = GROUP][order(-mean_pi)])
  fwrite(pw, file.path(TAB, "Data8_windowed_pi.csv"))
}

## ---------------------------------------------------------------------------
## F. Restorer versus maintainer lines.
## In sunflower hybrid breeding, R (restorer) and B (maintainer) lines are kept
## as separate pools, so they are expected to be differentiated. GRIN records
## the designation in the accession narrative, and USDA line names encode it
## (HA is a maintainer, RHA a restorer); scripts/grin_extract_type.py collects
## both. The designation is meaningful only for breeding lines.
## ---------------------------------------------------------------------------
if (!is.null(meta) && "HETEROTIC_GROUP" %in% names(meta)) {
  hg <- meta[HETEROTIC_GROUP %chin% c("R", "B") & SAMPLE %chin% indNames(gid_all)]
  cat(sprintf("\nheterotic group: %d R, %d B\n",
              hg[HETEROTIC_GROUP == "R", .N], hg[HETEROTIC_GROUP == "B", .N]))
  if (hg[, uniqueN(HETEROTIC_GROUP)] == 2 && hg[, min(table(HETEROTIC_GROUP))] >= 5) {
    g_rb  <- gid_all[indNames(gid_all) %chin% hg$SAMPLE]
    hf_rb <- genind2hierfstat(g_rb,
               pop = factor(hg[match(indNames(g_rb), SAMPLE), HETEROTIC_GROUP]))
    bs_rb <- basic.stats(hf_rb)
    print(bs_rb$overall)
    wc_rb <- pairwise.WCfst(hf_rb)
    print(wc_rb)
    ## strata must be set on this subset; the object inherits the cluster
    ## strata from the AMOVA above and would otherwise have no such level
    gc_rb <- as.genclone(g_rb)
    strata(gc_rb) <- data.frame(
      HETEROTIC_GROUP = factor(hg[match(indNames(gc_rb), SAMPLE), HETEROTIC_GROUP]))
    am_rb <- poppr.amova(gc_rb, ~HETEROTIC_GROUP, within = FALSE, quiet = TRUE)
    print(am_rb)
    set.seed(5); rb_sig <- randtest(am_rb, nrepet = 999); print(rb_sig)
    ov <- bs_rb$overall
    fwrite(rbind(
      data.table(MODEL = "Restorer versus maintainer",
                 SECTION = "Global F-statistics", TERM = names(ov),
                 DF = NA_real_, SUM_SQ = NA_real_, MEAN_SQ = NA_real_,
                 SIGMA = NA_real_, PERCENT = NA_real_,
                 STATISTIC = round(as.numeric(ov), 5), P_VALUE = NA_real_),
      data.table(MODEL = "Restorer versus maintainer",
                 SECTION = "Pairwise Weir and Cockerham F_ST",
                 TERM = "Maintainer versus restorer", DF = NA_real_,
                 SUM_SQ = NA_real_, MEAN_SQ = NA_real_, SIGMA = NA_real_,
                 PERCENT = NA_real_, STATISTIC = round(wc_rb[1, 2], 5),
                 P_VALUE = NA_real_),
      amova_tidy(am_rb, rb_sig, "Restorer versus maintainer")),
      file.path(TAB, "TableS12_restorer_vs_maintainer.csv"))
  } else {
    cat("too few labelled lines for the R versus B contrast\n")
  }
}

## ---------------------------------------------------------------------------
## G. Geographic origin, separated by the kind of origin event.
## The GRIN origin field mixes where an accession was collected, where it was
## bred, and who donated it. Only the collected subset carries geography in the
## sense a diversity study means, so REGION is tested there as well as overall.
## ---------------------------------------------------------------------------
if (!is.null(meta) && all(c("REGION", "ORIGIN_EVENT") %in% names(meta))) {
  ev <- meta[SAMPLE %chin% indNames(gid_all)]
  print(ev[, .N, by = .(ORIGIN_EVENT, REGION)][order(ORIGIN_EVENT, -N)])
  ev_rows <- list()
  for (e in c("Collected", "Developed")) {
    sub <- ev[ORIGIN_EVENT == e & REGION != "Unknown"]
    keep_r <- sub[, .N, by = REGION][N >= 5, REGION]
    sub <- sub[REGION %chin% keep_r]
    if (uniqueN(sub$REGION) < 2) next
    g_e <- gid_all[indNames(gid_all) %chin% sub$SAMPLE]
    hf_e <- genind2hierfstat(g_e, pop = factor(sub[match(indNames(g_e), SAMPLE), REGION]))
    cat(sprintf("\nREGION within ORIGIN_EVENT = %s: %d accessions, %d regions\n",
                e, nrow(sub), uniqueN(sub$REGION)))
    bs_e <- basic.stats(hf_e)$overall
    print(bs_e)
    ev_rows[[e]] <- data.table(ORIGIN_EVENT = e, N = nrow(sub),
                               N_REGIONS = uniqueN(sub$REGION),
                               as.data.table(as.list(round(bs_e, 4))))
  }
  if (length(ev_rows))
    fwrite(rbindlist(ev_rows), file.path(TAB, "TableS11c_origin_event_stats.csv"))
}

## Tajima's D is not computed here. The provider filtered both delivered SNP
## sets at MAF >= 5%, so the rare-allele tail that Watterson's theta depends on
## is absent from the data and D cannot be estimated. vcftools still writes
## grp_*.Tajima.D in step 05; those files are not used.

cat("step 14 done\n")
