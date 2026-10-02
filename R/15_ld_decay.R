## Step 15 - LD decay, Hill & Weir fit ------------------------------------------
source(file.path(Sys.getenv("PROJ"), "analysis/R/00_common.R"))

LD <- file.path(OUT, "05_ld")

## Reads either PopLDdecay *.stat(.gz) output (cols: Dist Mean_r^2 Mean_D' Sum_r^2 Sum_D' NumberPairs)
## or the plink fallback table written by slurm/06_ld_decay.slurm.
read_ld <- function(f) {
  if (grepl("\\.bin\\.tsv$", f)) {
    d <- fread(f)
    return(data.table(dist_bp = d$DIST_KB * 1000, r2 = d$MEAN_R2, n = d$N_PAIRS))
  }
  d <- fread(cmd = if (grepl("\\.gz$", f)) paste("zcat", shQuote(f)) else NULL,
             file = if (grepl("\\.gz$", f)) NULL else f)
  setnames(d, 1:2, c("Dist", "MeanR2"))
  np <- if (ncol(d) >= 6) d[[6]] else rep(1L, nrow(d))
  data.table(dist_bp = as.numeric(d$Dist), r2 = as.numeric(d$MeanR2), n = as.numeric(np))
}

files <- list.files(LD, pattern = "\\.stat(\\.gz)?$|\\.bin\\.tsv$", full.names = TRUE)
stopifnot(length(files) > 0)

ld <- rbindlist(lapply(files, function(f) {
  lab <- sub("\\.stat(\\.gz)?$|\\.bin\\.tsv$", "", basename(f))
  lab <- sub("^grp_", "", sub("^all$", "All accessions", lab))
  ## only the whole panel and the named clusters belong in the figure; any
  ## other stat file in the directory is a leftover and is skipped
  if (!grepl("^(All accessions|Q[0-9]+)$", lab)) return(NULL)
  d <- read_ld(f); d[, SET := lab]; d
}))
ld <- ld[dist_bp > 0 & is.finite(r2)]

## Hill & Weir (1988) E(r2). rho = C * d, C = 4*Ne*c per bp,
## n = number of sampled chromosomes.
hw <- function(C, d, n) {
  rho <- C * d
  ((10 + rho) / ((2 + rho) * (11 + rho))) *
    (1 + ((3 + rho) * (12 + 12 * rho + rho^2)) / (n * (2 + rho) * (11 + rho)))
}

nchr_sampled <- 2 * length(readLines(file.path(OUT, "01_qc", "keep_samples.txt")))

fits <- ld[, {
  st <- tryCatch(
    nls(r2 ~ hw(C, dist_bp, nchr_sampled), start = list(C = 3e-5),
        algorithm = "port", lower = c(C = 1e-9), upper = c(C = 1e-2),
        control = nls.control(maxiter = 500, warnOnly = TRUE)),
    error = function(e) NULL)
  if (is.null(st)) .(C = NA_real_) else .(C = unname(coef(st)["C"]))
}, by = SET]
print(fits)

## half-decay distance and the r2 = 0.2 / 0.1 thresholds, read off the fitted curve
grid <- data.table(dist_bp = seq(100, 2e6, by = 100))
curve <- fits[!is.na(C), {
  cc <- C
  copy(grid)[, r2 := hw(cc, dist_bp, nchr_sampled)]
}, by = SET]

thr <- curve[, .(
  fit_r2_at_1kb   = r2[which.min(abs(dist_bp - 1e3))],
  fit_r2_at_10kb  = r2[which.min(abs(dist_bp - 1e4))],
  fit_r2_at_100kb = r2[which.min(abs(dist_bp - 1e5))],
  fit_d_r2_0.2_kb = suppressWarnings(min(dist_bp[r2 <= 0.2], na.rm = TRUE) / 1000),
  fit_d_r2_0.1_kb = suppressWarnings(min(dist_bp[r2 <= 0.1], na.rm = TRUE) / 1000)
), by = SET]

## ---------------------------------------------------------------------------
## Empirical decay. The Hill-Weir model caps at E(r2) = 0.4545 at rho = 0 and
## decays smoothly, which cannot reproduce the shape of a tGBS marker set:
## SNPs cluster around restriction sites, so a large majority of pairs sit
## below 1 kb at high r2 and the curve then drops sharply. The fitted values
## above are kept for reference; the values reported in the manuscript are the
## binned empirical means and the crossings of a monotone fit to them.
## ---------------------------------------------------------------------------
BIN <- 2000
bin_curve <- function(x) {
  e <- x[, .(r2 = sum(r2 * n) / sum(n), pairs = sum(n)),
         by = .(SET, bin = floor(dist_bp / BIN))]
  e[, dist_kb := (bin + 0.5) * BIN / 1000]
  setorder(e, SET, bin)
  e[]
}

## isotonic (monotone non-increasing) fit, so a noisy bin cannot decide the
## crossing distance on its own. r2 at fixed distances is the stable statistic
## and is what the manuscript reports; the crossing distances are kept for
## comparison with the literature, but on a curve that plateaus a 0.03 change
## in r2 moves a crossing by tens of kb, so they are not quoted as primary
## results.
fixed_points <- function(e) e[, {
  o  <- order(bin)
  ir <- isoreg(bin[o], -r2[o])
  mono <- -ir$yf
  d <- (bin[o] + 0.5) * BIN / 1000
  f  <- function(t) { i <- which(mono <= t)[1]; if (is.na(i)) NA_real_ else d[i] }
  at <- function(x) mono[which.min(abs(d - x))]
  .(emp_r2_at_1kb   = at(1),   emp_r2_at_5kb   = at(5),
    emp_r2_at_10kb  = at(10),  emp_r2_at_20kb  = at(20),
    emp_r2_at_50kb  = at(50),  emp_r2_at_100kb = at(100),
    emp_r2_at_200kb = at(200), emp_r2_at_500kb = at(500),
    emp_d_r2_0.2_kb = f(0.20), emp_d_r2_0.15_kb = f(0.15),
    emp_d_r2_0.1_kb = f(0.10))
}, by = SET]

emp <- bin_curve(ld)
fwrite(emp[, .(SET, DISTANCE_KB = dist_kb, R2 = round(r2, 5), PAIRS = pairs)],
       file.path(TAB, "Data9_LD_binned_empirical.csv"))
cross <- fixed_points(emp)

## the Hill-Weir fit stays in the log as a reference; the manuscript table
## carries the empirical values only
tab6 <- Reduce(function(a, b) merge(a, b, by = "SET"), list(fits, thr, cross))
print(tab6)
tab6_out <- copy(cross)
setnames(tab6_out,
         c("emp_r2_at_1kb", "emp_r2_at_5kb", "emp_r2_at_10kb", "emp_r2_at_20kb",
           "emp_r2_at_50kb", "emp_r2_at_100kb", "emp_r2_at_200kb",
           "emp_r2_at_500kb", "emp_d_r2_0.2_kb", "emp_d_r2_0.15_kb",
           "emp_d_r2_0.1_kb"),
         c("R2_1KB", "R2_5KB", "R2_10KB", "R2_20KB", "R2_50KB", "R2_100KB",
           "R2_200KB", "R2_500KB", "DIST_R2_0.20_KB", "DIST_R2_0.15_KB",
           "DIST_R2_0.10_KB"))
num <- setdiff(names(tab6_out), "SET")
tab6_out[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(tab6_out, file.path(TAB, "Table3_LD_decay.csv"))

## ---------------------------------------------------------------------------
## Sensitivity of the whole-panel curve to imputation and to the heterozygosity
## filter. slurm/11_ld_sensitivity.slurm writes the sens_* files.
## ---------------------------------------------------------------------------
sfiles <- list.files(LD, pattern = "^sens_.*\\.stat(\\.gz)?$", full.names = TRUE)
if (length(sfiles)) {
  slab <- c(sens_nonimputed = "Non-imputed genotypes",
            sens_het020     = "Heterozygous call rate <= 0.20")
  lds <- rbindlist(lapply(sfiles, function(f) {
    key <- sub("\\.stat(\\.gz)?$", "", basename(f))
    d <- read_ld(f)
    d[, SET := if (!is.na(slab[key])) unname(slab[key]) else key]
    d
  }))
  lds <- lds[dist_bp > 0 & is.finite(r2)]
  cs <- fixed_points(bin_curve(lds))
  dcols <- grep("^emp_r2_at_", names(cross), value = TRUE)
  rows <- c(list(cross[SET == "All accessions"]), split(cs, seq_len(nrow(cs))))
  sens <- rbindlist(lapply(rows, function(r) data.table(
    SET = r$SET,
    DISTANCE_KB = as.numeric(sub("kb$", "", sub("^emp_r2_at_", "", dcols))),
    R2 = round(as.numeric(r[, ..dcols]), 4))))
  sens <- merge(sens,
                sens[SET == "All accessions", .(DISTANCE_KB, R2_REFERENCE = R2)],
                by = "DISTANCE_KB")
  sens[, `:=`(DIFFERENCE = round(R2 - R2_REFERENCE, 4),
              RATIO = round(R2 / R2_REFERENCE, 3))]
  setcolorder(sens, c("SET", "DISTANCE_KB", "R2", "R2_REFERENCE", "DIFFERENCE", "RATIO"))
  setorder(sens, SET, DISTANCE_KB)
  print(sens)
  fwrite(sens, file.path(TAB, "TableS23_LD_sensitivity.csv"))
}

## ---------------------------------------------------------------------------
## Figure
## ---------------------------------------------------------------------------
## The binned empirical means are plotted, not the Hill-Weir fit: the model
## cannot reproduce a curve that drops sharply and then plateaus.
p <- ggplot(emp[dist_kb <= 500], aes(dist_kb, r2, colour = SET)) +
  geom_line(linewidth = .7) +
  geom_hline(yintercept = 0.2, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = CLUSTER_COLS, name = NULL) +
  labs(x = "Physical distance (kb)",
       y = expression("Linkage disequilibrium," ~ r^2)) +
  theme_pub()
save_fig(p, "Fig9_ld_decay", 7.5, 5)

p2 <- p + scale_x_log10() + labs(x = "Physical distance (kb, log scale)")
save_fig(p2, "Fig9b_ld_decay_log", 7.5, 5)

## pair counts per bin, so the figure caption can state where the curve is
## supported by data
print(emp[, .(bins = .N, pairs = sum(pairs),
              pairs_under_5kb = sum(pairs[dist_kb < 5])), by = SET])


cat("step 15 done\n")
