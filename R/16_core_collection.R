## Step 16 - core collection sampling + kinship heat map -------------------------
## The smallest accession subset
## that retains the allelic richness and structure of the whole collection.
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
suppressPackageStartupMessages({ library(corehunter); library(rJava) })

PL  <- file.path(OUT, "01_qc")
CR  <- file.path(OUT, "07_core"); dir.create(CR, showWarnings = FALSE, recursive = TRUE)

## Core Hunter wants a biallelic 0/1/2 matrix.  Export it from the pruned set.
## (run once:  plink2 --bfile sunflower_pruned --export A --out sunflower_pruned_A )
raw <- file.path(PL, "sunflower_pruned_A.raw")
if (!file.exists(raw))
  stop("Run first:  plink2 --bfile ", PL, "/sunflower_pruned --allow-extra-chr --export A --out ",
       PL, "/sunflower_pruned_A")

d <- fread(raw)
ids <- d$IID
m <- as.matrix(d[, -(1:6)])
rownames(m) <- ids
m[is.na(m)] <- round(colMeans(m, na.rm = TRUE))[col(m)][is.na(m)]

gdata <- genotypes(data = m, format = "biparental")   # raw 0/1/2 dosages; taxa from rownames(m)

## Fraction of the collection to keep - report 10/15/20/25/30%
fractions <- c(0.10, 0.15, 0.20, 0.25, 0.30)
res <- rbindlist(lapply(fractions, function(fr) {
  set.seed(99)
  core <- sampleCore(gdata, size = fr,
                     obj = list(objective("EN", "MR"),      # entry-to-nearest-entry, modified Rogers
                                objective("SH")),           # Shannon diversity
                     mode = "default")
  writeLines(core$sel, file.path(CR, sprintf("core_%02d.txt", round(fr * 100))))
  ev <- evaluateCore(core$sel, gdata, objective("EN", "MR"))
  sh <- evaluateCore(core$sel, gdata, objective("SH"))
  he <- evaluateCore(core$sel, gdata, objective("HE"))
  data.table(FRACTION = fr, N = length(core$sel), EN_MR = ev, SHANNON = sh, HE = he)
}))

## whole-collection reference values
res <- rbind(res, data.table(FRACTION = 1, N = length(ids),
                             EN_MR   = evaluateCore(ids, gdata, objective("EN", "MR")),
                             SHANNON = evaluateCore(ids, gdata, objective("SH")),
                             HE      = evaluateCore(ids, gdata, objective("HE"))))
print(res)

asg <- fread(file.path(OUT, "03_admixture/group_assignment.txt"),
             header = FALSE, col.names = c("SAMPLE", "GROUP"))

## ---------------------------------------------------------------------------
## Stratified core. Core Hunter maximizes marker diversity over the whole panel
## and is free to drop an entire cluster, which makes the resulting subset
## unusable as a germplasm core. Allocating slots per cluster first, then
## optimizing within each cluster, keeps every group represented.
## ---------------------------------------------------------------------------
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
    cg <- sampleCore(gd, size = k,
                     obj = list(objective("EN", "MR"), objective("SH")),
                     mode = "default")
    sel <- c(sel, cg$sel)
  }
  sel
}

res_s <- rbindlist(lapply(fractions, function(fr) {
  sel <- strat_core(fr)
  writeLines(sel, file.path(CR, sprintf("core_stratified_%02d.txt", round(fr * 100))))
  data.table(FRACTION = fr, N = length(sel),
             EN_MR   = evaluateCore(sel, gdata, objective("EN", "MR")),
             SHANNON = evaluateCore(sel, gdata, objective("SH")),
             HE      = evaluateCore(sel, gdata, objective("HE")))
}))
res_s[, METHOD := "stratified"]
res[, METHOD := "Core Hunter, whole panel"]
tab7 <- rbind(res, res_s, fill = TRUE)[order(METHOD, FRACTION)]
fwrite(tab7, file.path(TAB, "Table3_core_collection.csv"))
print(tab7)

## cluster representation of both 20% cores
rep_tab <- function(sel, label) {
  merge(asg[, .(N_total = .N), by = GROUP],
        asg[SAMPLE %chin% sel, .(N_core = .N), by = GROUP],
        by = "GROUP", all.x = TRUE)[
          , .(GROUP, N_total, N_core = fifelse(is.na(N_core), 0L, N_core),
              PCT = round(100 * fifelse(is.na(N_core), 0L, N_core) / N_total, 1),
              METHOD = label)]
}
cov <- rbind(rep_tab(readLines(file.path(CR, "core_20.txt")), "Core Hunter, whole panel"),
             rep_tab(readLines(file.path(CR, "core_stratified_20.txt")), "stratified"))
print(dcast(cov, GROUP + N_total ~ METHOD, value.var = "PCT"))
fwrite(cov, file.path(TAB, "TableS15_core_cluster_representation.csv"))

p <- ggplot(tab7[FRACTION < 1], aes(100 * FRACTION, SHANNON, colour = METHOD)) +
  geom_line() + geom_point(size = 2.4) +
  geom_hline(yintercept = res[FRACTION == 1, SHANNON], linetype = 2, colour = "grey40") +
  scale_colour_manual(values = PAL[c(1, 3)], name = NULL,
                      labels = function(x) sub("^stratified$", "Stratified", x)) +
  labs(x = "Core size (% of collection)", y = "Shannon diversity index") + theme_pub()
save_fig(p, "Fig10_core_collection", 7, 4.5)

## ---------------------------------------------------------------------------
## Kinship heat map from the GCTA-form relationship matrix (plink2 --make-grm-list)
## ---------------------------------------------------------------------------
gf <- file.path(PL, "grm.grm")
if (file.exists(gf)) {
  gl <- fread(gf, col.names = c("i", "j", "nsnp", "k"))
  n  <- max(gl$i)
  K  <- matrix(0, n, n)
  K[cbind(gl$i, gl$j)] <- gl$k
  K[cbind(gl$j, gl$i)] <- gl$k
  gid_f <- fread(paste0(gf, ".id"), header = FALSE)     # plink2 writes grm.grm.id
  dimnames(K) <- list(gid_f$V2, gid_f$V2)
  ord <- hclust(as.dist(1 - K / max(K)), method = "average")$order
  kl <- as.data.table(reshape2::melt(K[ord, ord], varnames = c("A", "B"), value.name = "K"))
  p2 <- ggplot(kl, aes(A, B, fill = pmin(K, quantile(K, .995)))) + geom_raster() +
    scale_fill_viridis_c(option = "B", name = "Kinship") +
    labs(x = NULL, y = NULL) + theme_pub() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  save_fig(p2, "Fig11_kinship_heatmap", 7.5, 6.5)
  fwrite(as.data.table(K, keep.rownames = "SAMPLE"), file.path(TAB, "Data1_kinship.csv"))
}
cat("step 16 done\n")
