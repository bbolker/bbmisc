"""Fetch full texts for the sampled papers into papers/.

Order of sources: Europe PMC full-text XML (open access), OpenAlex open-access
PDF, then the McMaster EZproxy using the Firefox login session. Proxy requests
go one at a time with a pause between them, and the run stops on a login page
or repeated blocks. Every attempt is logged to fetch_log.csv; papers that
already have a file are skipped, so the script can be rerun.

Usage: python3 fetch_papers.py COOKIE_JAR CA_BUNDLE [DOI ...]
(DOIs default to every sampled paper in sample.csv)
"""

import csv
import html
import http.cookiejar
import re
import sys
import time
import urllib.parse
from datetime import datetime
from pathlib import Path

import requests

HERE = Path(__file__).parent
PAPERS = HERE / "papers"
LOG = HERE / "fetch_log.csv"
PROXY = "libaccess.lib.mcmaster.ca"
PAUSE = 20
UA = "Mozilla/5.0 (X11; Linux x86_64; rv:140.0) Gecko/20100101 Firefox/140.0"
BLOCK_TITLES = ("just a moment", "attention required", "access denied")


def safe_name(doi):
    return re.sub(r"[^A-Za-z0-9.]+", "_", doi)


def log(doi, source, status, detail=""):
    new = not LOG.exists()
    with LOG.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["time", "doi", "source", "status", "detail"])
        w.writerow([datetime.now().isoformat(timespec="seconds"), doi, source, status, detail])


def existing(doi):
    return [p for p in PAPERS.glob(safe_name(doi) + ".*") if p.suffix in (".pdf", ".xml")]


def page_title(text):
    m = re.search(r"<title[^>]*>(.*?)</title>", text, re.S | re.I)
    return html.unescape(m.group(1)).strip() if m else ""


class Blocked(Exception):
    pass


class LoggedOut(Exception):
    pass


def open_session():
    s = requests.Session()
    s.headers["User-Agent"] = "lit-survey-script/0.1"
    return s


def europepmc(s, doi):
    q = urllib.parse.quote(f'DOI:"{doi}"')
    r = s.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
              f"?query={q}&format=json&resultType=lite", timeout=60).json()
    res = r["resultList"]["result"]
    pmcid = res[0].get("pmcid") if res else None
    if not pmcid:
        return None
    x = s.get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML", timeout=120)
    if x.ok and len(x.content) > 20000 and b"<body" in x.content:
        return x.content, ".xml"
    return None


def openalex_pdf(s, doi):
    r = s.get(f"https://api.openalex.org/works/https://doi.org/{doi}"
              "?select=best_oa_location,open_access", timeout=60)
    if not r.ok:
        return None
    loc = r.json().get("best_oa_location") or {}
    ## best_oa_location can be a preprint (e.g. bioRxiv), whose analyses may
    ## differ from the published paper's
    if loc.get("version") != "publishedVersion":
        return None
    url = loc.get("pdf_url")
    if not url:
        return None
    p = s.get(url, timeout=120, headers={"User-Agent": UA})
    if p.ok and p.content[:5] == b"%PDF-":
        return p.content, ".pdf"
    return None


def proxied(url):
    """Rewrite a publisher URL to the proxy's host form."""
    u = urllib.parse.urlsplit(url)
    if u.netloc.endswith(PROXY):
        return url
    host = u.netloc.replace("-", "--").replace(".", "-") + "." + PROXY
    return urllib.parse.urlunsplit((u.scheme or "https", host, u.path, u.query, u.fragment))


class Proxy:
    def __init__(self, jar_path, ca_bundle):
        jar = http.cookiejar.MozillaCookieJar(jar_path)
        jar.load(ignore_discard=True, ignore_expires=True)
        for c in jar:
            c.expires = None  # expiry 0 means "session cookie" to curl, "expired" to Python
        self.s = requests.Session()
        self.s.cookies = jar
        self.s.headers["User-Agent"] = UA
        self.s.verify = ca_bundle
        self.n_blocks = 0
        self.last = 0.0

    def get(self, url, **kw):
        wait = PAUSE - (time.time() - self.last)
        if wait > 0:
            time.sleep(wait)
        r = self.s.get(url, timeout=120, **kw)
        self.last = time.time()
        title = page_title(r.text[:20000]) if r.content[:5] != b"%PDF-" else ""
        if "login to access library" in title.lower() or "/login" in urllib.parse.urlsplit(r.url).path:
            raise LoggedOut(url)
        if r.status_code in (403, 429) or title.lower().startswith(BLOCK_TITLES):
            self.n_blocks += 1
            if self.n_blocks >= 2:
                raise Blocked(f"{r.status_code} {title} at {r.url}")
        else:
            self.n_blocks = 0
        return r

    def fetch(self, doi):
        land = self.get(f"https://{PROXY}/login?url=https://doi.org/{doi}")
        if not urllib.parse.urlsplit(land.url).netloc.endswith(PROXY):
            land = self.get(proxied(land.url))
        if land.content[:5] == b"%PDF-":
            return land.content, ".pdf", land.url
        host = "https://" + urllib.parse.urlsplit(land.url).netloc
        candidates = []
        m = re.search(r'<meta[^>]+name="citation_pdf_url"[^>]+content="([^"]+)"', land.text)
        if m:
            candidates.append(proxied(html.unescape(m.group(1))))
        candidates += [f"{host}/doi/pdfdirect/{doi}?download=true", f"{host}/doi/pdf/{doi}"]
        for url in candidates:
            r = self.get(url, headers={"Referer": land.url})
            if r.content[:5] == b"%PDF-":
                return r.content, ".pdf", url
        return None, None, f"no PDF found from {land.url} (title: {page_title(land.text)[:60]})"


def main():
    jar, ca = sys.argv[1], sys.argv[2]
    dois = sys.argv[3:] or [r["doi"] for r in csv.DictReader(open(HERE / "sample.csv"))
                            if r["sampled"] == "y"]
    PAPERS.mkdir(exist_ok=True)
    oa = open_session()
    proxy = Proxy(jar, ca)
    for doi in dois:
        if existing(doi):
            continue
        got = None
        for name, fn in [("europepmc", europepmc), ("openalex", openalex_pdf)]:
            try:
                got = fn(oa, doi)
            except Exception as e:
                log(doi, name, "error", repr(e)[:200])
            if got:
                (PAPERS / (safe_name(doi) + got[1])).write_bytes(got[0])
                log(doi, name, "ok", f"{len(got[0])} bytes")
                print(f"{doi}: {name} ok", flush=True)
                break
        if got:
            continue
        if doi.startswith("10.1016/"):
            log(doi, "proxy", "skipped", "Elsevier: download manually")
            print(f"{doi}: Elsevier, needs manual download", flush=True)
            continue
        try:
            content, ext, where = proxy.fetch(doi)
        except LoggedOut:
            log(doi, "proxy", "logged out")
            print("Proxy session expired: log in again in Firefox and re-extract cookies.", flush=True)
            sys.exit(2)
        except Blocked as e:
            log(doi, "proxy", "blocked", str(e)[:200])
            print(f"Stopping: repeated blocks ({e})", flush=True)
            sys.exit(3)
        if content:
            (PAPERS / (safe_name(doi) + ext)).write_bytes(content)
            log(doi, "proxy", "ok", where)
            print(f"{doi}: proxy ok", flush=True)
        else:
            log(doi, "proxy", "failed", where)
            print(f"{doi}: proxy failed ({where})", flush=True)


if __name__ == "__main__":
    main()
