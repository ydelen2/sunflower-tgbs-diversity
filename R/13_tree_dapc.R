## Step 13 - NJ/BIONJ tree, bootstrap support, DAPC -----------------------------
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
suppressPackageStartupMessages({
  library(ape); library(adegenet); library(SNPRelate); library(ggtree)
})

PL  <- file.path(OUT, "01_qc")
TRD <- file.path(OUT, "06_tree"); dir.create(TRD, showWarnings = FALSE, recursive = TRUE)
gds <- file.path(TRD, "sunflower.gds")

## ---- GDS from the LD-pruned PLINK set -----------------------------------------
## Always rebuild: a cached GDS from an earlier marker set or accession list
## would be used silently and every downstream result would describe that set.
if (file.exists(gds)) unlink(gds)
snpgdsBED2GDS(paste0(PL, "/sunflower_pruned.bed"),
              paste0(PL, "/sunflower_pruned.fam"),
              paste0(PL, "/sunflower_pruned.bim"), gds, cvt.chr = "char")
g   <- snpgdsOpen(gds)
ids <- read.gdsn(index.gdsn(g, "sample.id"))
fam_n <- nrow(fread(paste0(PL, "/sunflower_pruned.fam"), header = FALSE))
stopifnot(length(ids) == fam_n)
cat(sprintf("tree and DAPC on %d accessions\n", length(ids)))

## ---- IBS distance and BIONJ ---------------------------------------------------
ibs <- snpgdsIBS(g, num.thread = 8, autosome.only = FALSE)
D   <- as.dist(1 - ibs$ibs)
attr(D, "Labels") <- ids
tr  <- bionj(D)
write.tree(tr, file.path(TRD, "sunflower_bionj.newick"))

## Bootstrap support over SNPs. Each replicate rebuilds a 274 x 49,523 distance
## matrix, so the number of replicates is set by BOOT_REPS (default 100) and can
## be lowered for a quick pass.
gmat <- snpgdsGetGeno(g, snp.id = ibs$snp.id, with.id = FALSE)
rownames(gmat) <- ids
set.seed(7)
bs <- boot.phylo(tr, gmat,
                 function(x) bionj(dist(x, method = "manhattan") / ncol(x)),
                 B = as.integer(Sys.getenv("BOOT_REPS", "100")),
                 trees = FALSE, quiet = TRUE)
tr$node.label <- bs
write.tree(tr, file.path(TRD, "sunflower_bionj_boot.newick"))
cat("median bootstrap support:", median(bs, na.rm = TRUE), "%\n")

## ---- plot, colored by ADMIXTURE cluster --------------------------------------
asg  <- fread(file.path(OUT, "03_admixture/group_assignment.txt"),
              header = FALSE, col.names = c("SAMPLE", "GROUP"))
meta <- read_meta()
if (!is.null(meta)) asg <- merge(asg, meta, by = "SAMPLE", all.x = TRUE)

p <- ggtree(tr, layout = "fan", open.angle = 12, linewidth = .3)
p <- p %<+% as.data.frame(asg) +
  geom_tippoint(aes(colour = GROUP), size = 1.6) +
  scale_colour_manual(values = CLUSTER_COLS, name = "Cluster") +
  theme(legend.position = "right")
save_fig(p, "Fig4_bionj_tree", 9, 9)

## rectangular version with bootstrap labels, for the supplement
p2 <- ggtree(tr, linewidth = .25) %<+% as.data.frame(asg) +
  geom_tippoint(aes(colour = GROUP), size = 1) +
  geom_tiplab(size = 1.1) +
  scale_colour_manual(values = CLUSTER_COLS, name = "Cluster")
save_fig(p2, "FigS6_bionj_tree_rect", 9, 26)

## ---- DAPC: a discriminant view, independent of ADMIXTURE ----------------------
gl <- new("genlight", gen = gmat, ind.names = ids, ploidy = 2)

set.seed(11)
grp <- find.clusters(gl, max.n.clust = 15, n.pca = 200,
                     choose.n.clust = FALSE, criterion = "diffNgroup")
nK  <- length(levels(grp$grp))

## then choose the number of PCs to retain by cross-validation rather than a guess
xv <- tryCatch(
  xvalDapc(tab(gl, NA.method = "mean"), grp = grp$grp,
           n.pca.max = 200, training.set = 0.9,
           result = "groupMean", n.rep = 10, parallel = "no"),
  error = function(e) { message("xvalDapc failed (", conditionMessage(e),
                                "); retaining 100 PCs"); NULL })
np <- if (is.null(xv)) 100L else as.integer(xv$DAPC$n.pca)
cat("find.clusters BIC-selected K =", nK, "; PCs retained =", np, "\n")

nda <- max(1, min(5, nK - 1))
dp <- dapc(gl, grp$grp, n.pca = np, n.da = nda)

## the settings selected, so that the numbers quoted in the text have a source
xs <- if (is.null(xv)) NULL else {
  ms <- xv$`Mean Successful Assignment by Number of PCs of PCA`
  rm_ <- xv$`Root Mean Squared Error by Number of PCs of PCA`
  k <- as.character(np)
  list(success = if (k %in% names(ms)) as.numeric(ms[[k]]) else NA_real_,
       rmse = if (k %in% names(rm_)) as.numeric(rm_[[k]]) else NA_real_)
}
settings <- data.table(
  SETTING = c("Number of discriminant groups (successive K-means, diffNgroup)",
              "Principal components available to K-means",
              "Principal components retained for the discriminant step",
              "Discriminant functions retained",
              "Cross-validation replicates", "Training fraction",
              "Mean proportion of correct assignment at the retained number of components",
              "Root mean squared error at the retained number of components"),
  VALUE = c(nK, 200, np, nda, 10, 0.9,
            if (is.null(xs)) NA_real_ else round(xs$success, 4),
            if (is.null(xs)) NA_real_ else round(xs$rmse, 4)))
print(settings)
fwrite(settings, file.path(TAB, "TableS2_dapc_settings.csv"))
dt <- data.table(SAMPLE = ids, DAPC = grp$grp,
                 LD1 = dp$ind.coord[, 1],
                 LD2 = if (ncol(dp$ind.coord) > 1) dp$ind.coord[, 2] else NA_real_)

## Points carry the ancestry-cluster color and ellipses outline the
## discriminant groups, so the agreement between the two methods is visible
## directly: each ellipse encloses points of a single color.
dt <- merge(dt, asg[, .(SAMPLE, GROUP)], by = "SAMPLE")
p3 <- ggplot(dt, aes(LD1, LD2)) +
  stat_ellipse(aes(group = DAPC), level = .67, linewidth = .4, colour = "grey40") +
  geom_point(aes(colour = GROUP), size = 2.2, alpha = .85) +
  scale_colour_manual(values = CLUSTER_COLS, name = "Cluster") +
  labs(x = "Discriminant function 1", y = "Discriminant function 2") + theme_pub()
save_fig(p3, "Fig5_dapc", 7, 5.5)

## BIC curve for the supplement
bic <- data.table(K = seq_along(grp$Kstat), BIC = as.numeric(grp$Kstat))
p4 <- ggplot(bic, aes(K, BIC)) + geom_line(colour = PAL[1]) +
  geom_point(colour = PAL[1]) + scale_x_continuous(breaks = bic$K) +
  labs(x = "Number of clusters (K)", y = "BIC") + theme_pub()
save_fig(p4, "FigS7_dapc_BIC", 6, 4)

fwrite(dt, file.path(TAB, "Data3_dapc_coordinates.csv"))

## concordance between DAPC and ADMIXTURE assignments
ct <- table(DAPC = dt$DAPC, ADMIXTURE = asg[match(dt$SAMPLE, asg$SAMPLE), GROUP])
print(ct)
fwrite(as.data.table(as.data.frame.matrix(ct), keep.rownames = "DAPC_GROUP"),
       file.path(TAB, "TableS7_dapc_vs_admixture.csv"))

snpgdsClose(g)
cat("step 13 done\n")
