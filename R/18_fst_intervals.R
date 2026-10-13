## Step 18 - uncertainty of the pairwise F_ST estimates ---------------------------
## Per-locus Weir and Cockerham variance components for every pair of ancestry
## clusters, on the full non-imputed filtered marker set, then a
## delete-one-chromosome jackknife (which respects linkage) and a bootstrap over
## loci. The point estimate is the ratio of summed components, the same
## quantity VCFtools reports as the weighted F_ST in Table S10.
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))
suppressPackageStartupMessages({ library(vcfR); library(hierfstat) })

asg <- fread(file.path(OUT, "03_admixture/group_assignment_pops.txt"),
             header = FALSE, col.names = c("SAMPLE", "GROUP"))
t3  <- fread(file.path(TAB, "pairwise_fst_clusters.csv"))

vcf <- read.vcfR(file.path(OUT, "00_vcf/sunflower.filtered.vcf.gz"), verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT")
chr <- vcf@fix[, "CHROM"]
rm(vcf)
## genotype strings to hierfstat two-digit codes; unphased and phased forms
code <- c("0/0" = 11L, "0/1" = 12L, "1/0" = 12L, "1/1" = 22L,
          "0|0" = 11L, "0|1" = 12L, "1|0" = 12L, "1|1" = 22L)
G <- matrix(code[gt], nrow = nrow(gt), dimnames = dimnames(gt))
rm(gt)
G <- t(G)                       # individuals in rows, loci in columns
stopifnot(all(asg$SAMPLE %in% rownames(G)))

## per-locus variance components (a: among populations, b: among individuals
## within populations, w: within individuals), in chunks of loci
components <- function(ids, pop, chunk = 20000L) {
  L <- ncol(G)
  out <- vector("list", ceiling(L / chunk))
  for (i in seq_along(out)) {
    cols <- ((i - 1L) * chunk + 1L):min(i * chunk, L)
    dat <- data.frame(pop = factor(pop), G[ids, cols, drop = FALSE])
    w <- wc(dat, diploid = TRUE)
    sl <- as.data.frame(w$sigma.loc)
    ia <- grep("a$", names(sl)); ib <- grep("b$", names(sl)); iw <- grep("w$", names(sl))
    if (length(ia) != 1L || length(ib) != 1L || length(iw) != 1L) { ia <- 1L; ib <- 2L; iw <- 3L }
    out[[i]] <- data.table(A = sl[[ia]], B = sl[[ib]], W = sl[[iw]])
  }
  rbindlist(out)
}

set.seed(18)
NBOOT <- 1000L
groups <- sort(unique(asg$GROUP))
pairs  <- combn(groups, 2, simplify = FALSE)
res <- rbindlist(lapply(pairs, function(pr) {
  ids <- asg[GROUP %in% pr, SAMPLE]
  pop <- asg[match(ids, SAMPLE), GROUP]
  cmp <- components(ids, pop)
  cmp[, TOT := A + B + W]
  ok  <- cmp[, is.finite(A) & is.finite(TOT)]
  a <- cmp$A[ok]; tt <- cmp$TOT[ok]; ch <- chr[ok]
  fst <- sum(a) / sum(tt)
  ## delete-one-chromosome jackknife
  chs <- sort(unique(ch)); m <- length(chs)
  theta <- vapply(chs, function(cc) { k <- ch != cc; sum(a[k]) / sum(tt[k]) }, numeric(1))
  se_jack <- sqrt((m - 1) / m * sum((theta - mean(theta))^2))
  ## bootstrap over loci
  L <- length(a)
  boot <- vapply(seq_len(NBOOT), function(b) { k <- sample.int(L, L, replace = TRUE); sum(a[k]) / sum(tt[k]) }, numeric(1))
  v <- t3[(POP1 == pr[1] & POP2 == pr[2]) | (POP1 == pr[2] & POP2 == pr[1])]
  cat(sprintf("%s-%s  F_ST %.4f (VCFtools %.4f)  jackknife SE %.4f  bootstrap CI %.4f-%.4f\n",
              pr[1], pr[2], fst, v$WEIGHTED_FST[1], se_jack, quantile(boot, .025), quantile(boot, .975)))
  data.table(POP1 = pr[1], POP2 = pr[2],
             N1 = sum(pop == pr[1]), N2 = sum(pop == pr[2]), N_LOCI = L,
             FST = round(fst, 4),
             FST_VCFTOOLS_WEIGHTED = round(v$WEIGHTED_FST[1], 4),
             JACKKNIFE_SE = round(se_jack, 5),
             CI_LOW_JACKKNIFE = round(fst - 1.96 * se_jack, 4),
             CI_HIGH_JACKKNIFE = round(fst + 1.96 * se_jack, 4),
             CI_LOW_BOOTSTRAP = round(quantile(boot, .025), 4),
             CI_HIGH_BOOTSTRAP = round(quantile(boot, .975), 4),
             N_CHROMOSOMES = m, N_BOOTSTRAP = NBOOT)
}))
print(res)
fwrite(res, file.path(TAB, "TableS10b_fst_intervals.csv"))
cat("step 18 done\n")
