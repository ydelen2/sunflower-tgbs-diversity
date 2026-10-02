## Shared setup for all R steps -------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

PROJ <- Sys.getenv("PROJ", unset = file.path(Sys.getenv("WORK"), "sunflower_tgbs"))
OUT  <- file.path(PROJ, "results")
ANA  <- file.path(PROJ, "analysis")
FIG  <- file.path(OUT, "08_figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
TAB  <- file.path(OUT, "09_tables");  dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

CHR_LAB <- setNames(paste0("Ha", 1:17), as.character(1:17))

theme_pub <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = .25, colour = "grey90"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          legend.key = element_blank())
}

## 12-color qualitative palette, color-blind safe enough for K<=12
PAL <- c("#4C72B0","#DD8452","#55A868","#C44E52","#8172B3","#937860",
         "#DA8BC3","#8C8C8C","#CCB974","#64B5CD","#B07AA1","#76B7B2")

## One palette for the ancestry clusters everywhere. Positional palettes shift
## when "admixed" sorts first as a factor level, which put Q1 in a different
## color in the tree than in the ancestry barplot.
CLUSTER_COLS <- c(Q1 = PAL[1], Q2 = PAL[2], Q3 = PAL[3], Q4 = PAL[4],
                  Q5 = PAL[5], admixed = "grey60",
                  `All accessions` = "grey25")

save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 400)
  invisible(p)
}

read_meta <- function() {
  f <- file.path(ANA, "metadata", "sample_metadata.csv")
  if (!file.exists(f)) {
    warning("sample_metadata.csv not found - a priori grouping analyses will be skipped")
    return(NULL)
  }
  ## Every column is read as character. Several label columns are blank for the
  ## first hundred rows, and type guessing can turn the later values into NA.
  m <- fread(f, colClasses = "character", na.strings = c("", "NA"))
  stopifnot("SAMPLE" %in% names(m))
  for (cl in intersect(c("OIL_PCT", "YEAR", "LAT", "LON"), names(m)))
    m[, (cl) := suppressWarnings(as.numeric(get(cl)))]
  m
}
