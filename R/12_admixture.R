## Step 12 - ADMIXTURE: CV curve, Q barplots, cluster assignment ----------------
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))

AD  <- file.path(OUT, "03_admixture")
fam <- fread(file.path(OUT, "01_qc/sunflower_pruned.fam"), header = FALSE)
ids <- fam$V2

cv <- fread(file.path(AD, "cv_summary.tsv"))
cv <- cv[!is.na(CV_ERROR)]

## --- CV error across the K grid ------------------------------------------------
agg <- cv[, .(mean = mean(CV_ERROR), sd = sd(CV_ERROR), min = min(CV_ERROR)), by = K][order(K)]
bestK <- agg[which.min(mean), K]
cat("Lowest mean CV error at K =", bestK, "\n")

p1 <- ggplot(agg, aes(K, mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .2, colour = "grey60") +
  geom_line(colour = PAL[1]) + geom_point(size = 2.4, colour = PAL[1]) +
  geom_point(data = agg[K == bestK], colour = PAL[4], size = 3.6) +
  scale_x_continuous(breaks = agg$K) +
  labs(x = "Number of ancestral clusters (K)",
       y = sprintf("Cross-validation error (mean \u00b1 SD, %d runs)", max(cv$REP))) +
  theme_pub()
save_fig(p1, "Fig3a_admixture_cv", 6.5, 4)
fwrite(agg, file.path(TAB, "TableS5_admixture_cv.csv"))

## --- Q matrices ----------------------------------------------------------------
best <- fread(file.path(AD, "best_rep_per_K.tsv"), header = FALSE,
              col.names = c("K", "REP", "CV"))
readQ <- function(K, rep) {
  f <- file.path(AD, sprintf("K%d_rep%d", K, rep), sprintf("sunflower_pruned.%d.Q", K))
  q <- as.matrix(fread(f)); rownames(q) <- ids; q
}

meta <- read_meta()
qlist <- list()
for (i in seq_len(nrow(best))) {
  K <- best$K[i]; if (K < 2) next
  qlist[[as.character(K)]] <- readQ(K, best$REP[i])
}

## order accessions once, by the best-K solution, so every barplot is comparable
Qb  <- qlist[[as.character(bestK)]]
asg <- data.table(SAMPLE = ids,
                  CLUSTER = paste0("Q", max.col(Qb)),
                  MAXQ    = apply(Qb, 1, max))
asg[, GROUP := fifelse(MAXQ >= 0.60, CLUSTER, "admixed")]
ord <- order(asg$CLUSTER, -asg$MAXQ)

long <- rbindlist(lapply(names(qlist), function(k) {
  q <- qlist[[k]]
  data.table(SAMPLE = rep(ids, ncol(q)),
             K = as.integer(k),
             CL = rep(paste0("Q", seq_len(ncol(q))), each = length(ids)),
             Q = as.vector(q))
}))
long[, SAMPLE := factor(SAMPLE, levels = ids[ord])]

p2 <- ggplot(long, aes(SAMPLE, Q, fill = CL)) +
  geom_col(width = 1) +
  facet_grid(K ~ ., switch = "y") +
  scale_fill_manual(values = rep(PAL, 3), guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = sprintf("%d accessions", length(ids)), y = "Ancestry proportion") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing.y = unit(1, "pt"), strip.placement = "outside")
save_fig(p2, "Fig3b_admixture_barplots_K2toKmax", 11, 1 + 0.9 * length(qlist))

## the single best-K panel, at publication size
p3 <- ggplot(long[K == bestK], aes(SAMPLE, Q, fill = CL)) +
  geom_col(width = 1) +
  scale_fill_manual(values = PAL, name = "Cluster") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry proportion", title = sprintf("K = %d", bestK)) +
  theme_pub() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p3, sprintf("Fig3c_admixture_K%d", bestK), 11, 3)

## --- group assignment consumed by steps 05 (Fst) and 06 (LD) -------------------
fwrite(asg[, .(SAMPLE, GROUP)], file.path(AD, "group_assignment.txt"),
       sep = "\t", col.names = FALSE)

## Clusters only, for pairwise F_ST / per-group diversity / per-group LD.
fwrite(asg[GROUP != "admixed", .(SAMPLE, GROUP)],
       file.path(AD, "group_assignment_pops.txt"), sep = "\t", col.names = FALSE)
cat(sprintf("assigned to a cluster: %d | admixed: %d\n",
            asg[GROUP != "admixed", .N], asg[GROUP == "admixed", .N]))

## K = 3 assignment, for comparison with the k-means three subpopulations
## reported in the published GWAS papers on this panel.
if ("3" %in% names(qlist)) {
  Q3  <- qlist[["3"]]
  a3  <- data.table(SAMPLE = ids, K3 = paste0("G", max.col(Q3)),
                    MAXQ3 = apply(Q3, 1, max))
  a3[, K3_GROUP := fifelse(MAXQ3 >= 0.60, K3, "admixed")]
  ct3 <- table(K3 = a3$K3_GROUP, bestK = asg$GROUP)
  cat("\nK=3 versus K=", bestK, " assignment:\n", sep = "")
  print(ct3)
  fwrite(as.data.table(as.data.frame.matrix(ct3), keep.rownames = "K3_GROUP"),
         file.path(TAB, "TableS7b_K3_vs_bestK.csv"))
  fwrite(a3[, .(SAMPLE, K3_GROUP, MAXQ3)], file.path(TAB, "Data6_K3_assignment.csv"))
}

Qout <- as.data.table(Qb)
setnames(Qout, paste0("Q", seq_len(ncol(Qout))))
fwrite(cbind(asg, Qout), file.path(TAB, "TableS6_admixture_Q.csv"))
print(asg[, .N, by = GROUP][order(-N)])
cat("step 12 done\n")
