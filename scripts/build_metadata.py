#!/usr/bin/env python3
"""
Python equivalent of R/01_build_metadata.R - joins the GRIN scrape onto the
sample list and writes metadata/sample_metadata.csv.  Runs anywhere Python does
(the R version exists for the Swan workflow; both produce the same file).

  python3 scripts/build_metadata.py
"""
import csv, os, re, sys

MD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "metadata")

REGION = {
    "North America": {"UNITED STATES", "CANADA", "MEXICO"},
    "South America": {"ARGENTINA", "CHILE", "PERU", "BOLIVIA", "URUGUAY", "BRAZIL"},
    "Eastern Europe / W Asia": {"RUSSIAN FEDERATION", "UKRAINE", "ROMANIA", "BULGARIA",
                                "HUNGARY", "SERBIA", "MOLDOVA", "POLAND", "TURKEY",
                                "FORMER SOVIET UNION", "YUGOSLAVIA", "CZECHOSLOVAKIA",
                                "CZECH REPUBLIC", "SLOVAKIA", "GEORGIA", "ARMENIA",
                                "AZERBAIJAN", "KAZAKHSTAN", "UZBEKISTAN", "SOVIET UNION"},
    "Western Europe": {"FRANCE", "SPAIN", "ITALY", "GERMANY", "NETHERLANDS", "PORTUGAL",
                       "GREECE", "AUSTRIA", "SWITZERLAND", "UNITED KINGDOM", "BELGIUM",
                       "DENMARK", "SWEDEN", "NORWAY", "FINLAND"},
    "Asia / Oceania": {"CHINA", "INDIA", "JAPAN", "PAKISTAN", "IRAN", "AUSTRALIA",
                       "ISRAEL", "SYRIA", "IRAQ", "AFGHANISTAN", "MONGOLIA", "KOREA, SOUTH",
                       "KOREA, NORTH", "THAILAND", "MYANMAR"},
    "Africa": {"SOUTH AFRICA", "EGYPT", "MOROCCO", "KENYA", "ZIMBABWE", "TANZANIA",
               "ETHIOPIA", "ALGERIA", "TUNISIA", "SUDAN"},
}


def country_of(raw: str) -> str:
    """GRIN origin strings are often 'State/Province, Country' - keep the country."""
    c = (raw or "").strip()
    if c.upper().startswith("FORMER, SOVIET"):
        return "Former Soviet Union"
    if "," in c:
        c = c.rsplit(",", 1)[1].strip()
    return c


def region_of(country: str) -> str:
    c = country_of(country).upper()
    if not c or c in ("UNKNOWN", "UNCERTAIN"):
        return "Unknown"
    if "SOVIET" in c:
        return "Eastern Europe / W Asia"
    if "SERBIA" in c or "MONTENEGRO" in c:
        return "Eastern Europe / W Asia"
    if c == "JORDAN":
        return "Asia / Oceania"
    for r, s in REGION.items():
        if c in s:
            return r
    return "Other"


def status_of(raw: str) -> str:
    s = (raw or "").lower()
    if "wild" in s:
        return "Wild"
    if "landrace" in s:
        return "Landrace"
    if "genetic" in s:
        return "Genetic stock"
    if any(k in s for k in ("breed", "line", "inbred")):
        return "Breeding line"
    if "cultivat" in s or "cultivar" in s:
        return "Cultivar"
    return "Unknown"


def main() -> int:
    samples = [l.strip() for l in open(os.path.join(MD, "samples_287.txt")) if l.strip()]

    grin = {}
    pf = os.path.join(MD, "grin_passport.csv")
    if os.path.exists(pf):
        for row in csv.DictReader(open(pf, encoding="utf-8")):
            grin[row["ACCESSION"]] = row
    else:
        print("grin_passport.csv missing - run scripts/grin_fetch.py first", file=sys.stderr)

    ## optional oil-content classifier from the observation table
    oil, seed_white = {}, {}
    of = os.path.join(MD, "grin_observations.csv")
    if os.path.exists(of):
        for row in csv.DictReader(open(of, encoding="utf-8")):
            acc = row["ACCESSION"]
            if row["TRAIT"] == "PER_OIL":
                m = re.match(r"[\d.]+", row["VALUE"])
                if m:
                    oil.setdefault(acc, []).append(float(m.group(0)))
            if row["TRAIT"] in ("SEED_COLOR", "SEED_COLOR1"):
                if re.search(r"WHITE|GREY|GRAY", row["VALUE"]):
                    seed_white[acc] = True

    ## preserve hand-edits from a previous run
    prev = {}
    outfile = os.path.join(MD, "sample_metadata.csv")
    if os.path.exists(outfile):
        for row in csv.DictReader(open(outfile, encoding="utf-8")):
            prev[row["SAMPLE"]] = row

    cols = ["SAMPLE", "ACCESSION", "SUBSAMPLE", "NAME", "TAXON", "COUNTRY", "COUNTRY_RAW", "REGION",
            "STATUS", "STATUS_RAW", "ORIGIN_EVENT", "TYPE", "TYPE_AUTO", "OIL_PCT",
            "HETEROTIC_GROUP", "YEAR", "SOURCE_HISTORY"]
    n_country = n_status = 0
    rows = []
    for s in samples:
        num = re.sub(r"^PI_", "", s)
        sub = ""
        m = re.search(r"(C\d+)$", num)
        if m:
            sub = m.group(1)
            num = num[: -len(sub)]
        acc = f"PI {num}"
        g = grin.get(acc, {})
        o = oil.get(acc)
        mean_oil = sum(o) / len(o) if o else None
        t_auto = ""
        if mean_oil is not None:
            if mean_oil >= 40:
                t_auto = "oilseed"
            elif mean_oil < 32 and seed_white.get(acc):
                t_auto = "confection"
        p = prev.get(s, {})
        row = {
            "SAMPLE": s, "ACCESSION": acc, "SUBSAMPLE": sub,
            "NAME": g.get("NAME", ""), "TAXON": g.get("TAXON", ""),
            "COUNTRY": country_of(g.get("COUNTRY", "")),
            "COUNTRY_RAW": g.get("COUNTRY", ""),
            "REGION": region_of(g.get("COUNTRY", "")),
            "STATUS": status_of(g.get("STATUS_RAW", "")), "STATUS_RAW": g.get("STATUS_RAW", ""),
            "ORIGIN_EVENT": g.get("ORIGIN_EVENT", ""),
            "TYPE": p.get("TYPE", ""),                       # hand-curated, preserved
            "TYPE_AUTO": t_auto,
            "OIL_PCT": f"{mean_oil:.1f}" if mean_oil is not None else "",
            "HETEROTIC_GROUP": p.get("HETEROTIC_GROUP", ""), # hand-curated, preserved
            "YEAR": g.get("YEAR", ""), "SOURCE_HISTORY": g.get("SOURCE_HISTORY", ""),
        }
        if row["COUNTRY"]:
            n_country += 1
        if row["STATUS"] != "Unknown":
            n_status += 1
        rows.append(row)

    with open(outfile, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    from collections import Counter
    print(f"wrote {outfile}  ({len(rows)} samples)", file=sys.stderr)
    print(f"  with country : {n_country}", file=sys.stderr)
    print(f"  with status  : {n_status}", file=sys.stderr)
    print("  REGION       :", dict(Counter(r['REGION'] for r in rows)), file=sys.stderr)
    print("  STATUS       :", dict(Counter(r['STATUS'] for r in rows)), file=sys.stderr)
    print("  TYPE_AUTO    :", dict(Counter(r['TYPE_AUTO'] or 'unclassified' for r in rows)),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
