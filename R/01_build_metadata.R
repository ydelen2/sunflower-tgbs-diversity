## Step 01 - build sample_metadata.csv from the GRIN-Global scrape ---------------
##
## Input: metadata/grin_passport.csv, produced by scripts/grin_fetch.py, which
## pulls each accession's public GRIN-Global detail page (passport + phenotype
## observations).  Fallback: a manually downloaded metadata/grin_download.csv
## from https://npgsweb.ars-grin.gov/gringlobal/search (the query list is in
## metadata/grin_accession_query.txt).
##
## Anything GRIN does not have (oil vs confection type, maintainer heterotic
## group R/B) must be filled by hand in the generated file; the script keeps
## hand-edits on re-run.

suppressPackageStartupMessages(library(data.table))

ANA <- file.path(Sys.getenv("PROJ"), "analysis")
MD  <- file.path(ANA, "metadata")

samples <- readLines(file.path(MD, "samples_287.txt"))
core    <- data.table(SAMPLE = samples)
core[, PI := sub("^PI_", "", SAMPLE)]
core[, ACCESSION := paste("PI", sub("C[0-9]+$", "", PI))]   # PI_490282C2 -> "PI 490282"
core[, SUBSAMPLE := ifelse(grepl("C[0-9]+$", PI), sub(".*(C[0-9]+)$", "\\1", PI), "")]

writeLines(unique(core$ACCESSION), file.path(MD, "grin_accession_query.txt"))

pf <- file.path(MD, "grin_passport.csv")
gf <- file.path(MD, "grin_download.csv")
if (file.exists(pf)) {
  g <- fread(pf)
  g <- g[, .(ACCESSION, TAXON, NAME, COUNTRY, ORIGIN_EVENT, STATUS = STATUS_RAW,
             LIFE_CYCLE, YEAR, SOURCE_HISTORY)]
} else if (file.exists(gf)) {
  g <- fread(gf)
  ## GRIN column names vary between download templates; normalize the common ones
  ren <- c("Accession"="ACCESSION", "Accession Number"="ACCESSION", "accession_number"="ACCESSION",
           "Origin"="COUNTRY", "Country"="COUNTRY", "origin"="COUNTRY",
           "State"="STATE", "Improvement Status"="STATUS", "improvement_status"="STATUS",
           "Plant Name"="NAME", "Name"="NAME",
           "Latitude"="LAT", "Longitude"="LON",
           "Received Date"="YEAR", "Taxonomy"="TAXON", "Taxon"="TAXON")
  for (k in intersect(names(ren), names(g))) setnames(g, k, ren[[k]])
  keep <- intersect(c("ACCESSION","COUNTRY","STATE","STATUS","NAME","LAT","LON","YEAR","TAXON"), names(g))
  g <- unique(g[, ..keep])
} else {
  message("Neither grin_passport.csv nor grin_download.csv found.\n",
          "Run:  python3 scripts/grin_fetch.py --samples metadata/samples_287.txt ",
          "--cache metadata/grin_cache --out metadata")
  tmpl <- copy(core)
  tmpl[, `:=`(COUNTRY = NA_character_, REGION = NA_character_, STATE = NA_character_,
              STATUS = NA_character_, TYPE = NA_character_, HETEROTIC_GROUP = NA_character_,
              NAME = NA_character_, LAT = NA_real_, LON = NA_real_, YEAR = NA_integer_)]
  fwrite(tmpl, file.path(MD, "sample_metadata.csv"))
  quit(save = "no")
}

out <- merge(core, g, by = "ACCESSION", all.x = TRUE)

## ---- geographic region, for AMOVA / Mantel ------------------------------------
region_map <- function(x) {
  x <- toupper(trimws(as.character(x)))
  fcase(
    x %chin% c("UNITED STATES","CANADA","MEXICO"),                       "North America",
    x %chin% c("ARGENTINA","CHILE","PERU","BOLIVIA","URUGUAY","BRAZIL"), "South America",
    x %chin% c("RUSSIAN FEDERATION","UKRAINE","ROMANIA","BULGARIA","HUNGARY","SERBIA",
               "MOLDOVA","POLAND","TURKEY","FORMER SOVIET UNION","YUGOSLAVIA"), "Eastern Europe / W Asia",
    x %chin% c("FRANCE","SPAIN","ITALY","GERMANY","NETHERLANDS","PORTUGAL","GREECE"), "Western Europe",
    x %chin% c("CHINA","INDIA","JAPAN","PAKISTAN","IRAN","AUSTRALIA"),   "Asia / Oceania",
    x %chin% c("SOUTH AFRICA","EGYPT","MOROCCO","KENYA","ZIMBABWE"),     "Africa",
    is.na(x) | x == "" | x == "UNKNOWN",                                  "Unknown",
    default = "Other")
}
out[, REGION := region_map(COUNTRY)]

## ---- improvement status, collapsed to three levels ---------------------------
## GRIN's phrasing on the public pages: "Wild material", "Landrace",
## "Cultivated material", "Cultivar", "Breeding material", "Landrace",
## "Genetic material" (genetic stocks, kept separate from breeding lines)
out[, STATUS := fcase(
  grepl("wild",      STATUS, ignore.case = TRUE), "Wild",
  grepl("landrace",  STATUS, ignore.case = TRUE), "Landrace",
  grepl("genetic",             STATUS, ignore.case = TRUE), "Genetic stock",
  grepl("breed|line|inbred",   STATUS, ignore.case = TRUE), "Breeding line",
  grepl("cultivar|cultivated", STATUS, ignore.case = TRUE), "Cultivar",
  default = "Unknown")]

## Market type and heterotic group are not structured GRIN fields, but curators
## record them in the accession narrative and USDA line names encode the group
## (HA maintainer, RHA restorer). scripts/grin_extract_type.py collects both
## from the cached pages; run it before this script.
tgf <- file.path(MD, "grin_type_group.csv")
if (file.exists(tgf)) {
  tg <- fread(tgf)[, .(ACCESSION, TYPE, TYPE_SOURCE, HETEROTIC_GROUP, GROUP_SOURCE)]
  tg[tg == ""] <- NA
  out <- merge(out, tg, by = "ACCESSION", all.x = TRUE)
} else {
  for (cl in c("TYPE", "TYPE_SOURCE", "HETEROTIC_GROUP", "GROUP_SOURCE"))
    out[, (cl) := NA_character_]
}

## preserve manual edits from a previous run
prev <- file.path(MD, "sample_metadata.csv")
if (file.exists(prev)) {
  old <- fread(prev)
  for (cl in intersect(c("TYPE", "HETEROTIC_GROUP"), names(old))) {
    out[old, on = "SAMPLE", (cl) := fifelse(is.na(get(cl)), get(paste0("i.", cl)), get(cl))]
  }
}

setcolorder(out, c("SAMPLE","ACCESSION","SUBSAMPLE","NAME","COUNTRY","REGION","STATE",
                   "STATUS","TYPE","HETEROTIC_GROUP","LAT","LON","YEAR"))
fwrite(out, prev)

cat("wrote", prev, "\n")
print(out[, .N, by = REGION][order(-N)])
print(out[, .N, by = STATUS][order(-N)])
cat("\nStill missing COUNTRY:", out[is.na(COUNTRY) | COUNTRY == "", .N], "accessions\n")
cat("Still missing TYPE    :", out[is.na(TYPE), .N], "accessions - fill by hand\n")
