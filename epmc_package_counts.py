import requests
import time
import csv

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# package name -> list of surface-form variants to OR together
PACKAGES = {
    "DESeq2": ['"DESeq2"'],
    "edgeR": ['"edgeR"'],
    "ANCOM-BC": ['"ANCOM-BC"', '"ANCOMBC"', '"ANCOM BC"'],
    "ALDEx2": ['"ALDEx2"'],
    "corncob": ['"corncob"'],
    "NBZIMM": ['"NBZIMM"'],
    "Maaslin2": ['"MaAsLin2"', '"Maaslin2"'],
    "Maaslin3": ['"MaAsLin3"', '"Maaslin3"'],
    "gllvm": ['"gllvm"'],
    "mvabund": ['"mvabund"'],
}

MICROBIOME_TERMS = ['"16S"', '"16S rRNA"', "amplicon", "microbiome"]

YEARS = range(2015, 2026)

FIELD_TO_QUERY_FIELDS = {
    "fulltext_count": (None, None),
    "mixed_count": (None, "ABSTRACT"),
    "abstract_count": ("ABSTRACT", "ABSTRACT"),
}


def hit_count(query, tries=5):
    # Europe PMC occasionally returns a well-formed response with a spuriously
    # low hitCount under rapid successive requests (observed even with rate
    # limiting). Retrying with backoff resolves it; the caller additionally
    # cross-checks the three tiers for monotonicity (see get_consistent_triplet).
    params = {
        "query": query,
        "format": "json",
        "pageSize": 1,
    }
    for attempt in range(tries):
        try:
            r = requests.get(BASE, params=params, timeout=30)
            r.raise_for_status()
            data = r.json()
            if "hitCount" in data:
                return data["hitCount"]
        except Exception as e:
            print(f"  error (attempt {attempt + 1}): {query[:80]}... -> {e}")
        time.sleep(1.5 * (attempt + 1))
    print(f"  FAILED after {tries} tries: {query[:80]}...")
    return None


def build_query(pkg_variants, year, pkg_field, micro_field):
    pkg_clause = " OR ".join(
        f"{pkg_field}:{v}" if pkg_field else v for v in pkg_variants
    )
    micro_clause = " OR ".join(
        f"{micro_field}:{t}" if micro_field else t for t in MICROBIOME_TERMS
    )
    return f"({pkg_clause}) AND ({micro_clause}) AND PUB_YEAR:{year}"


def get_consistent_triplet(pkg, variants, year, max_rounds=6):
    # mixed_count is a subset of fulltext_count's match set, and abstract_count
    # is a subset of mixed_count's, so a valid triplet must satisfy
    # abstract <= mixed <= fulltext. A transient API hiccup on any one of the
    # three calls breaks that ordering, so re-issue all three until consistent.
    ft = mixed = ab = None
    for round_ in range(max_rounds):
        ft = hit_count(build_query(variants, year, pkg_field=None, micro_field=None))
        time.sleep(0.5)
        mixed = hit_count(build_query(variants, year, pkg_field=None, micro_field="ABSTRACT"))
        time.sleep(0.5)
        ab = hit_count(build_query(variants, year, pkg_field="ABSTRACT", micro_field="ABSTRACT"))
        time.sleep(0.5)
        if None not in (ft, mixed, ab) and ab <= mixed <= ft:
            return ft, mixed, ab
        print(f"  retry {pkg} {year}: got ft={ft} mixed={mixed} ab={ab} (inconsistent), round {round_ + 1}")
        time.sleep(2)
    print(f"  GAVE UP on {pkg} {year}, keeping last values ft={ft} mixed={mixed} ab={ab}")
    return ft, mixed, ab


def reverify_zero(pkg, variants, year, field, n_independent=3):
    # A 0 can be a stable-but-wrong response (observed: identical across
    # several immediate retries within one run, yet resolved cleanly by a
    # fresh query issued later) rather than a genuine absence, so a 0 gets
    # re-queried independently rather than trusted on the first pass.
    pkg_field, micro_field = FIELD_TO_QUERY_FIELDS[field]
    query = build_query(variants, year, pkg_field, micro_field)
    results = [hit_count(query) for _ in range(n_independent)]
    for _ in range(n_independent - 1):
        time.sleep(0.6)
    clean = [v for v in results if v is not None]
    return max(clean) if clean else 0


rows = []
for pkg, variants in PACKAGES.items():
    for year in YEARS:
        triplet = dict(zip(
            ("fulltext_count", "mixed_count", "abstract_count"),
            get_consistent_triplet(pkg, variants, year),
        ))

        for field in ("fulltext_count", "mixed_count", "abstract_count"):
            if triplet[field] == 0:
                best = reverify_zero(pkg, variants, year, field)
                if best != 0:
                    print(f"  {pkg} {year} {field}: 0 -> {best} on independent re-check")
                    triplet[field] = best

        # re-enforce monotonicity after any zero was bumped
        triplet["mixed_count"] = min(triplet["mixed_count"], triplet["fulltext_count"])
        triplet["abstract_count"] = min(triplet["abstract_count"], triplet["mixed_count"])

        print(f"{pkg}\t{year}\tfulltext={triplet['fulltext_count']}\tmixed={triplet['mixed_count']}\tabstract={triplet['abstract_count']}")
        rows.append({"package": pkg, "year": year, **triplet})

out_path = "/home/bolker/Documents/students/agronah/package_counts_by_year.csv"
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["package", "year", "fulltext_count", "mixed_count", "abstract_count"])
    writer.writeheader()
    writer.writerows(rows)

print(f"\nWrote {len(rows)} rows to {out_path}")
