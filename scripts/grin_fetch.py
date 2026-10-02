#!/usr/bin/env python3
"""
Pull passport data AND phenotype observations for a list of USDA PI accessions
from the GRIN-Global public website.

The public site resolves an accession directly from its identifier, so
"PI_162675" -> "PI 162675" -> one page, one request:

    https://npgsweb.ars-grin.gov/gringlobal/accessiondetail?accid=PI%20162675

Every page is cached under --cache, so re-runs cost nothing and the USDA server
is hit once per accession.  Default rate is 1 request/second; do not lower it.

Usage
  python3 grin_fetch.py --samples ../metadata/samples_287.txt \
                        --cache  ../metadata/grin_cache \
                        --out    ../metadata
Outputs
  grin_passport.csv      one row per accession
  grin_observations.csv  long format: ACCESSION, CATEGORY, TRAIT, DESCRIPTION, VALUE
  grin_failed.txt        accessions whose page could not be parsed
"""
import argparse, csv, html, os, re, sys, time, urllib.parse, urllib.request

UA = "Mozilla/5.0 (compatible; academic germplasm study; contact via USDA GRIN feedback)"
BASE = "https://npgsweb.ars-grin.gov/gringlobal/accessiondetail?accid="

# stable ASP.NET label ids on the accession detail page
FIELDS = {
    "TAXON":      "MainContent_ctrlSum_lblTaxon",
    "NAME":       "MainContent_ctrlSum_lblTopName",
    "ORIGIN_RAW": "MainContent_ctrlSum_lblOrigin",
    "SITE":       "MainContent_ctrlSum_lblSite",
    "RECEIVED":   "MainContent_ctrlSum_lblReceived",
    "PI_ASSIGNED":"MainContent_ctrlSum_lblAssigned",
    "STATUS_RAW": "MainContent_ctrlSum_lblStatus",
    "LIFE_FORM":  "MainContent_ctrlSum_lblLifeForm",
    "LIFE_CYCLE": "MainContent_ctrlSum_lblLifeCycle",
    "BACKUP":     "MainContent_ctrlSum_lblBackup",
    "GENUS":      "MainContent_ctrlTaxon_lblGenus",
    "FAMILY":     "MainContent_ctrlTaxon_lblFamily",
}


def clean(s: str) -> str:
    s = re.sub(r"<[^>]+>", " ", s)
    s = html.unescape(s)
    s = s.replace("\ufffd", "-")
    return re.sub(r"\s+", " ", s).strip()


def fetch(acc: str, cache: str, delay: float, retries: int = 3) -> str | None:
    path = os.path.join(cache, acc.replace(" ", "_") + ".html")
    if os.path.exists(path) and os.path.getsize(path) > 5000:
        return open(path, encoding="utf-8", errors="replace").read()
    url = BASE + urllib.parse.quote(acc)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            body = urllib.request.urlopen(req, timeout=45).read().decode("utf-8", "replace")
            open(path, "w", encoding="utf-8").write(body)
            time.sleep(delay)
            return body
        except Exception as e:                       # noqa: BLE001
            print(f"  [{acc}] attempt {attempt+1}/{retries}: {type(e).__name__} {e}",
                  file=sys.stderr)
            time.sleep(delay * (attempt + 2))
    return None


def parse_passport(acc: str, h: str) -> dict | None:
    row = {"ACCESSION": acc}
    got = 0
    for key, sid in FIELDS.items():
        m = re.search(r'<span id="%s"[^>]*>(.*?)</span>' % re.escape(sid), h, re.S)
        row[key] = clean(m.group(1)) if m else ""
        if row[key]:
            got += 1
    if got < 3:
        return None

    # "Collected - Argentina" / "Developed - United States" / "Unknown"
    o = row.get("ORIGIN_RAW", "")
    if "-" in o:
        head, _, tail = o.partition("-")
        row["ORIGIN_EVENT"] = head.strip()          # Collected / Developed / Bred / Donated
        row["COUNTRY"] = tail.strip()
    else:
        row["ORIGIN_EVENT"] = ""
        row["COUNTRY"] = o.strip()

    # received date -> year
    m = re.search(r"(19|20)\d{2}", row.get("RECEIVED", ""))
    row["YEAR"] = m.group(0) if m else ""

    # source-history block often names the collecting locality
    m = re.search(r"Source History(.*?)Accession Names", h, re.S)
    row["SOURCE_HISTORY"] = clean(m.group(1))[:500] if m else ""
    return row


OBS_RE = re.compile(
    r"<tr>\s*<td>([A-Z ]+)</td>\s*"
    r"<td><a href=\"descriptordetail\.aspx\?id=\d+\"[^>]*>([A-Z0-9_]+)</a></td>\s*"
    r"<td>([^<]*)</td>\s*<td>([^<]*)</td>")


def parse_observations(acc: str, h: str) -> list[dict]:
    out = []
    for cat, trait, desc, val in OBS_RE.findall(h):
        v = clean(val)
        if v:
            out.append({"ACCESSION": acc, "CATEGORY": cat.strip(),
                        "TRAIT": trait, "DESCRIPTION": clean(desc), "VALUE": v})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True,
                    help="file with sample IDs, one per line (PI_162675 or 'PI 162675')")
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--delay", type=float, default=1.0, help="seconds between requests")
    a = ap.parse_args()

    os.makedirs(a.cache, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)

    # PI_490282C2 and PI_490282C3 are re-samples of accession PI 490282
    accs, seen = [], set()
    for ln in open(a.samples):
        s = ln.strip()
        if not s:
            continue
        num = re.sub(r"^PI[_ ]?", "", s)
        num = re.sub(r"C\d+$", "", num)             # drop the C2/C3 subsample suffix
        acc = f"PI {num}"
        if acc not in seen:
            seen.add(acc); accs.append(acc)

    print(f"{len(accs)} unique accessions to fetch (delay {a.delay}s)", file=sys.stderr)

    passport, obs, failed = [], [], []
    for i, acc in enumerate(accs, 1):
        h = fetch(acc, a.cache, a.delay)
        if h is None:
            failed.append(acc); continue
        p = parse_passport(acc, h)
        if p is None:
            failed.append(acc); continue
        passport.append(p)
        obs.extend(parse_observations(acc, h))
        if i % 25 == 0:
            print(f"  {i}/{len(accs)}  ({len(obs)} observations so far)", file=sys.stderr)

    pcols = ["ACCESSION", "TAXON", "GENUS", "FAMILY", "NAME", "COUNTRY", "ORIGIN_EVENT",
             "ORIGIN_RAW", "STATUS_RAW", "LIFE_FORM", "LIFE_CYCLE", "SITE", "BACKUP",
             "RECEIVED", "YEAR", "PI_ASSIGNED", "SOURCE_HISTORY"]
    with open(os.path.join(a.out, "grin_passport.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=pcols, extrasaction="ignore")
        w.writeheader(); w.writerows(passport)

    with open(os.path.join(a.out, "grin_observations.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["ACCESSION", "CATEGORY", "TRAIT", "DESCRIPTION", "VALUE"])
        w.writeheader(); w.writerows(obs)

    if failed:
        open(os.path.join(a.out, "grin_failed.txt"), "w").write("\n".join(failed) + "\n")

    print(f"\npassport rows     : {len(passport)}", file=sys.stderr)
    print(f"observation rows  : {len(obs)}", file=sys.stderr)
    print(f"failed            : {len(failed)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
