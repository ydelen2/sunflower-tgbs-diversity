#!/usr/bin/env python3
"""
Extract market type (oilseed / confection) and heterotic group (restorer /
maintainer) from the cached GRIN-Global accession pages.

Neither is a structured GRIN field, but curators record both in the narrative
and description free text, and USDA sunflower line names encode the heterotic
group: HA is a maintainer (B) line, RHA is a restorer (R) line.

In GRIN wording "non-oilseed" and "nonoil" denote the confection type.

Usage
  python3 grin_extract_type.py --cache ../metadata/grin_cache \
                               --out   ../metadata/grin_type_group.csv
"""
import argparse, csv, html, pathlib, re, sys

SPANS = ("MainContent_ctrlNarrative_lblNarrative",)
BLOCK = re.compile(r'class="(?:panel-body|col-md-10)"[^>]*>(.*?)</div>', re.S)


def text_of(page: str) -> str:
    parts = []
    for sid in SPANS:
        m = re.search(r'<span id="%s"[^>]*>(.*?)</span>' % re.escape(sid), page, re.S)
        if m:
            parts.append(m.group(1))
    parts.extend(BLOCK.findall(page))
    t = " ".join(parts)
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", html.unescape(t)).strip()


def name_of(page: str) -> str:
    m = re.search(r'<span id="MainContent_ctrlSum_lblTopName"[^>]*>(.*?)</span>', page, re.S)
    if not m:
        return ""
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html.unescape(m.group(1)))).strip()


RE_NONOIL = re.compile(r"non-?oil(seed)?", re.I)
RE_OIL = re.compile(r"\boil-?seed\b", re.I)
RE_REST = re.compile(r"\brestorer\b", re.I)
RE_MAINT = re.compile(r"\bmaintainer\b", re.I)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = []
    for f in sorted(pathlib.Path(a.cache).glob("*.html")):
        page = f.read_text(encoding="utf-8", errors="replace")
        acc = f.stem.replace("_", " ")
        txt = text_of(page)
        nm = name_of(page)

        if RE_NONOIL.search(txt):
            typ, typ_src = "confection", "narrative"
        elif RE_OIL.search(txt):
            typ, typ_src = "oilseed", "narrative"
        else:
            typ, typ_src = "", ""

        if RE_REST.search(txt):
            grp, grp_src = "R", "narrative"
        elif RE_MAINT.search(txt):
            grp, grp_src = "B", "narrative"
        elif re.match(r"^RHA[\s-]*\d", nm, re.I):
            grp, grp_src = "R", "line name"
        elif re.match(r"^HA[\s-]*\d", nm, re.I):
            grp, grp_src = "B", "line name"
        else:
            grp, grp_src = "", ""

        rows.append({"ACCESSION": acc, "NAME": nm,
                     "TYPE": typ, "TYPE_SOURCE": typ_src,
                     "HETEROTIC_GROUP": grp, "GROUP_SOURCE": grp_src})

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["ACCESSION", "NAME", "TYPE", "TYPE_SOURCE",
                                           "HETEROTIC_GROUP", "GROUP_SOURCE"])
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    for col in ("TYPE", "HETEROTIC_GROUP"):
        counts = {}
        for r in rows:
            counts[r[col] or "(none)"] = counts.get(r[col] or "(none)", 0) + 1
        print(f"{col}: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())),
              file=sys.stderr)
    print(f"{n} accession pages parsed, wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
