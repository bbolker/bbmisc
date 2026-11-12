"""Build toc_sept2026.csv: every research article in the September 2026 issues of
the five journals surveyed by Forstmeier & Schielzeth (2011), with abstracts.

Article lists come from Crossref (print-issue date in September 2026) for the
four print journals. Proc R Soc B is online-only and not yet indexed in Crossref
for September, so its September issues (by issue number) come from OpenAlex.
Abstracts come from Crossref when present, otherwise OpenAlex, otherwise
Europe PMC; missing ones are left blank and filled in by hand.
"""

import csv
import html
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

UA = {"User-Agent": "lit-survey-script/0.1"}
OUT = Path(__file__).with_name("toc_sept2026.csv")

PRINT_JOURNALS = {
    "AmNat": "0003-0147",
    "AnimBehav": "0003-3472",
    "Ecology": "0012-9658",
    "Evolution": "0014-3820",
}
PROCB_ISSN = "0962-8452"
PROCB_ISSUES = ["2078", "2079"]


def get_json(url):
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.load(urllib.request.urlopen(req, timeout=120))
        except Exception:
            if attempt == 2:
                raise
            time.sleep(5)


def clean_jats(s):
    """Strip JATS/HTML tags from a Crossref abstract."""
    s = re.sub(r"<[^>]+>", " ", s or "")
    return re.sub(r"\s+", " ", html.unescape(s)).strip()


def openalex_abstract(inv):
    """Rebuild abstract text from OpenAlex's inverted index."""
    if not inv:
        return ""
    pos = {i: w for w, idx in inv.items() for i in idx}
    return " ".join(pos[i] for i in sorted(pos))


def crossref_rows(journal, issn):
    flt = ("from-print-pub-date:2026-09-01,until-print-pub-date:2026-09-30,"
           "type:journal-article")
    url = (f"https://api.crossref.org/journals/{issn}/works?filter={flt}&rows=200"
           "&select=DOI,title,volume,issue,page,abstract,subject")
    rows = []
    for it in get_json(url)["message"]["items"]:
        rows.append({
            "journal": journal,
            "doi": it["DOI"].lower(),
            "volume": it.get("volume", ""),
            "issue": it.get("issue", ""),
            "page": it.get("page", ""),
            "title": clean_jats((it.get("title") or [""])[0]),
            "abstract": clean_jats(it.get("abstract", "")),
            "abstract_source": "crossref" if it.get("abstract") else "",
        })
    return rows


def procb_rows():
    flt = (f"primary_location.source.issn:{PROCB_ISSN},"
           f"biblio.issue:{'|'.join(PROCB_ISSUES)}")
    url = ("https://api.openalex.org/works?per-page=200"
           "&select=doi,title,biblio,abstract_inverted_index,type"
           f"&filter={flt}")
    rows = []
    for r in get_json(url)["results"]:
        if not r.get("doi"):
            continue
        b = r["biblio"]
        ab = openalex_abstract(r.get("abstract_inverted_index"))
        rows.append({
            "journal": "ProcB",
            "doi": r["doi"].replace("https://doi.org/", "").lower(),
            "volume": b.get("volume") or "",
            "issue": b.get("issue") or "",
            "page": b.get("first_page") or "",
            "title": clean_jats(r.get("title") or ""),
            "abstract": ab,
            "abstract_source": "openalex" if ab else "",
        })
    return rows


def fill_openalex(rows):
    todo = [r for r in rows if not r["abstract"]]
    for i in range(0, len(todo), 40):
        chunk = todo[i:i + 40]
        q = "|".join("https://doi.org/" + r["doi"] for r in chunk)
        url = ("https://api.openalex.org/works?per-page=50"
               "&select=doi,abstract_inverted_index&filter=doi:"
               + urllib.parse.quote(q, safe=":/|"))
        found = {res["doi"].replace("https://doi.org/", "").lower():
                 openalex_abstract(res.get("abstract_inverted_index"))
                 for res in get_json(url)["results"]}
        for r in chunk:
            if found.get(r["doi"]):
                r["abstract"], r["abstract_source"] = found[r["doi"]], "openalex"


def fill_europepmc(rows):
    for r in rows:
        if r["abstract"]:
            continue
        q = urllib.parse.quote(f'DOI:"{r["doi"]}"')
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
               f"query={q}&format=json&resultType=core")
        res = get_json(url)["resultList"]["result"]
        if res and res[0].get("abstractText"):
            r["abstract"] = clean_jats(res[0]["abstractText"])
            r["abstract_source"] = "europepmc"
        time.sleep(0.5)


def main():
    rows = []
    for journal, issn in PRINT_JOURNALS.items():
        rows += crossref_rows(journal, issn)
    rows += procb_rows()
    fill_openalex(rows)
    fill_europepmc(rows)

    fields = ["journal", "doi", "volume", "issue", "page", "title",
              "abstract_source", "abstract"]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    for j in [*PRINT_JOURNALS, "ProcB"]:
        sub = [r for r in rows if r["journal"] == j]
        print(f"{j:10s} {len(sub):3d} articles, "
              f"{sum(bool(r['abstract']) for r in sub):3d} with abstracts")


if __name__ == "__main__":
    main()
