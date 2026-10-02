## Step 19 - LD decay per chromosome and phenotypic coverage of the cores ----------
## (a) whole-panel LD decay per chromosome, from slurm/12_ld_by_chromosome.slurm,
##     with marker density and chromosome length from Table 1 and the correlations
##     of r-squared at 100 kb with both
## (b) phenotypic coverage of the two 20 percent cores against GRIN descriptors
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))

## ---- (a) LD decay by chromosome ------------------------------------------------
LDC <- file.path(OUT, "05_ld", "bychr")
files <- list.files(LDC, pattern = "^chr_Ha[0-9]+\\.stat(\\.gz)?$", full.names = TRUE)

read_ld <- function(f) {
  d <- fread(cmd = if (grepl("\\.gz$", f)) paste("zcat", shQuote(f)) else NULL,
             file = if (grepl("\\.gz$", f)) NULL else f)
  setnames(d, 1:2, c("Dist", "MeanR2"))
  np <- if (ncol(d) >= 6) d[[6]] else rep(1L, nrow(d))
  data.table(dist_bp = as.numeric(d$Dist), r2 = as.numeric(d$MeanR2), n = as.numeric(np))
}

BIN <- 2000
bin_curve <- function(x) {
  e <- x[, .(r2 = sum(r2 * n) / sum(n), pairs = sum(n)),
         by = .(SET, bin = floor(dist_bp / BIN))]
  setorder(e, SET, bin)
  e[]
}
fixed_points <- function(e) e[, {
  o  <- order(bin)
  ir <- isoreg(bin[o], -r2[o])
  mono <- -ir$yf
  d <- (bin[o] + 0.5) * BIN / 1000
  f  <- function(t) { i <- which(mono <= t)[1]; if (is.na(i)) NA_real_ else d[i] }
  at <- function(x) mono[which.min(abs(d - x))]
  .(PAIRS = sum(pairs),
    R2_1KB = at(1), R2_5KB = at(5), R2_10KB = at(10), R2_20KB = at(20),
    R2_50KB = at(50), R2_100KB = at(100), R2_200KB = at(200), R2_500KB = at(500),
    DIST_R2_0.20_KB = f(0.20), DIST_R2_0.15_KB = f(0.15), DIST_R2_0.10_KB = f(0.10))
}, by = SET]

if (length(files)) {
  ld <- rbindlist(lapply(files, function(f) {
    d <- read_ld(f)
    d[, SET := sub("\\.stat(\\.gz)?$", "", sub("^chr_", "", basename(f)))]
    d
  }))
  ld <- ld[dist_bp > 0 & is.finite(r2)]
  tab <- fixed_points(bin_curve(ld))
  setnames(tab, "SET", "CHROMOSOME")
  ## marker density and length from the marker summary table
  t1 <- fread(file.path(TAB, "Table1_marker_summary.csv"))
  t1 <- t1[CHRL %in% tab$CHROMOSOME, .(CHROMOSOME = CHRL, LENGTH_MB = round(LEN_Mb, 1),
                                        SNPS = N_SNP, SNPS_PER_MB = round(SNP_per_Mb, 1))]
  tab <- merge(t1, tab, by = "CHROMOSOME")
  tab[, ORD := as.integer(sub("Ha", "", CHROMOSOME))]
  setorder(tab, ORD); tab[, ORD := NULL]
  num <- setdiff(names(tab), c("CHROMOSOME", "PAIRS", "SNPS", "LENGTH_MB", "SNPS_PER_MB"))
  tab[, (num) := lapply(.SD, round, 3), .SDcols = num]
  print(tab)
  fwrite(tab, file.path(TAB, "TableS18_LD_by_chromosome.csv"))

  ## is the between-chromosome variation explained by marker density or length?
  ct <- function(x, y, lab) {
    r <- cor.test(x, y)
    data.table(COMPARISON = lab, PEARSON_R = round(unname(r$estimate), 3),
               N_CHROMOSOMES = length(x), P_VALUE = signif(r$p.value, 3))
  }
  cors <- rbind(
    ct(tab$SNPS_PER_MB, tab$R2_100KB, "r-squared at 100 kb against SNPs per Mb"),
    ct(tab$LENGTH_MB, tab$R2_100KB, "r-squared at 100 kb against chromosome length"),
    ct(tab$SNPS_PER_MB, tab$R2_10KB, "r-squared at 10 kb against SNPs per Mb"),
    ct(tab$LENGTH_MB, tab$R2_10KB, "r-squared at 10 kb against chromosome length"))
  print(cors)
  fwrite(cors, file.path(TAB, "TableS18b_LD_by_chromosome_correlations.csv"))
} else {
  message("no per-chromosome stat files under ", LDC, "; run slurm/12 first")
}

## ---- (b) phenotypic coverage of the cores ------------------------------------
CR  <- file.path(OUT, "07_core")
acc <- fread(file.path(TAB, "TableS1_accession_list.csv"), colClasses = "character")
obs <- fread(file.path(ANA, "metadata", "grin_observations.csv"), colClasses = "character")
obs <- obs[ACCESSION %in% acc$ACCESSION]

cores <- list(
  `Whole panel` = readLines(file.path(CR, "core_20.txt")),
  `Stratified`  = readLines(file.path(CR, "core_stratified_20.txt")))
cores <- lapply(cores, function(s) acc[SAMPLE %in% s, ACCESSION])
stopifnot(all(lengths(cores) == 55L))

continuous <- c(PER_OIL = "Seed oil content (per cent)",
                OLEIC = "Oleic acid (per cent of oil)",
                LINOLEIC = "Linoleic acid (per cent of oil)",
                DAY_FLOW = "Days to flowering",
                BRANCHING1 = "Plants without branching (per cent)",
                SEED_SIZE1 = "Seed passing a 20/64 inch screen (per cent)",
                SCLER_BASALSTEM_ROT = "Plants resistant to Sclerotinia basal stem rot (per cent)",
                PHOMOPSIS = "Plants resistant to Phomopsis (per cent)")
categorical <- c(SEED_COLOR = "Seed colour", SEED_STRIPES = "Seed stripes",
                 SEED_MOTTLING = "Seed mottling", SEED_PIG_ANTH = "Seed anthocyanin",
                 SEED_STRIPE_COLOR = "Seed stripe colour", SEED_THICK = "Seed thickness",
                 RAY_FLOW_COLOR = "Ray flower colour", RAY_FLOW_NUM = "Ray flower number",
                 POLLEN_COLOR = "Pollen colour", STIGMA_COLOR = "Stigma colour",
                 BRACT_TIGHT = "Bract tightness", DAY_FLOW_UNIF = "Flowering uniformity")

## one value per accession: mean of numeric observations; modal category code
## (the leading integer of the GRIN code, mixture percentages in brackets dropped)
num_val <- function(v) suppressWarnings(as.numeric(sub("\\s*\\(.*\\)$", "", v)))
cat_val <- function(v) {
  code <- sub("^\\s*(\\d+).*$", "\\1", sub("\\s*\\(\\d+%\\)$", "", v))
  code[!grepl("^\\d+$", code)] <- NA
  code
}
mode1 <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else names(sort(table(x), decreasing = TRUE))[1] }

per_acc_num <- obs[TRAIT %in% names(continuous), .(VALUE = mean(num_val(VALUE), na.rm = TRUE)),
                   by = .(TRAIT, ACCESSION)][is.finite(VALUE)]
per_acc_cat <- obs[TRAIT %in% names(categorical), .(VALUE = mode1(cat_val(VALUE))),
                   by = .(TRAIT, ACCESSION)][!is.na(VALUE)]

rows <- list()
for (tr in names(continuous)) {
  pan <- per_acc_num[TRAIT == tr]
  if (nrow(pan) < 30) next
  q <- quantile(pan$VALUE, probs = seq(0, 1, 0.1), na.rm = TRUE)
  for (sc in names(cores)) {
    cv <- pan[ACCESSION %in% cores[[sc]], VALUE]
    dec <- if (length(cv)) length(unique(pmin(findInterval(cv, q, rightmost.closed = TRUE), 10))) else 0L
    ks <- if (length(cv) >= 5) suppressWarnings(ks.test(cv, pan$VALUE)$p.value) else NA_real_
    rows[[length(rows) + 1L]] <- data.table(
      DESCRIPTOR = tr, DESCRIPTION = continuous[[tr]], TYPE = "continuous", SCHEME = sc,
      PANEL_N = nrow(pan), CORE_N = length(cv),
      PANEL_MIN = min(pan$VALUE), PANEL_MAX = max(pan$VALUE),
      PANEL_MEAN = round(mean(pan$VALUE), 2), CORE_MEAN = if (length(cv)) round(mean(cv), 2) else NA_real_,
      CORE_MIN = if (length(cv)) min(cv) else NA_real_, CORE_MAX = if (length(cv)) max(cv) else NA_real_,
      RANGE_COVERAGE = if (length(cv)) round((max(cv) - min(cv)) / (max(pan$VALUE) - min(pan$VALUE)), 3) else NA_real_,
      DECILES_COVERED = dec, KS_P = signif(ks, 3),
      PANEL_CATEGORIES = NA_integer_, CORE_CATEGORIES = NA_integer_, CATEGORY_COVERAGE = NA_real_,
      CATEGORIES_MISSED = NA_character_, MISSED_CATEGORY_PANEL_N = NA_character_)
  }
}
for (tr in names(categorical)) {
  pan <- per_acc_cat[TRAIT == tr]
  if (nrow(pan) < 30) next
  pc <- sort(unique(pan$VALUE))
  freq <- table(pan$VALUE)
  for (sc in names(cores)) {
    cv <- pan[ACCESSION %in% cores[[sc]], VALUE]
    cc <- sort(unique(cv))
    miss <- setdiff(pc, cc)
    rows[[length(rows) + 1L]] <- data.table(
      DESCRIPTOR = tr, DESCRIPTION = categorical[[tr]], TYPE = "categorical", SCHEME = sc,
      PANEL_N = nrow(pan), CORE_N = length(cv),
      PANEL_MIN = NA_real_, PANEL_MAX = NA_real_, PANEL_MEAN = NA_real_, CORE_MEAN = NA_real_,
      CORE_MIN = NA_real_, CORE_MAX = NA_real_,
      RANGE_COVERAGE = NA_real_, DECILES_COVERED = NA_integer_, KS_P = NA_real_,
      PANEL_CATEGORIES = length(pc), CORE_CATEGORIES = length(cc),
      CATEGORY_COVERAGE = round(length(cc) / length(pc), 3),
      CATEGORIES_MISSED = paste(miss, collapse = ";"),
      MISSED_CATEGORY_PANEL_N = paste(as.integer(freq[miss]), collapse = ";"))
  }
}
cov <- rbindlist(rows)
print(cov[, .(DESCRIPTOR, SCHEME, PANEL_N, CORE_N, PANEL_MEAN, CORE_MEAN, RANGE_COVERAGE,
              DECILES_COVERED, KS_P, CATEGORY_COVERAGE, MISSED_CATEGORY_PANEL_N)])
cat("\nsummary by scheme\n")
print(cov[, .(mean_range_coverage = round(mean(RANGE_COVERAGE, na.rm = TRUE), 3),
              min_range_coverage = round(min(RANGE_COVERAGE, na.rm = TRUE), 3),
              mean_deciles = round(mean(DECILES_COVERED, na.rm = TRUE), 2),
              min_deciles = min(DECILES_COVERED, na.rm = TRUE),
              ks_below_0.05 = sum(KS_P < 0.05, na.rm = TRUE),
              mean_category_coverage = round(mean(CATEGORY_COVERAGE, na.rm = TRUE), 3),
              min_category_coverage = round(min(CATEGORY_COVERAGE, na.rm = TRUE), 3)), by = SCHEME])
fwrite(cov, file.path(TAB, "TableS20_core_phenotype_coverage.csv"))
cat("step 19 done\n")
