"""Convert a downloaded paper (PDF or Europe PMC XML) to text and print the
passages relevant to coding: the statistical-methods section plus sentences
matching model/selection/correction keywords. Reference-list lines are dropped.

Usage: python3 snippets.py papers/<file>.pdf|.xml [max_chars]
"""

import html
import re
import subprocess
import sys
from pathlib import Path

KEYWORDS = re.compile(
    r"(?i)\b(statistic\w*|model\w*|glmm?s?|lmms?|mixed[- ]effects?|lme4|glmmTMB|nlme|brms|"
    r"MCMCglmm|bayes\w*|gams?|mgcv|interaction\w*|AICc?|BIC|stepwise|backward|forward|"
    r"drop1|dredge|MuMIn|model selection|model averag\w*|simplif\w*|removed|non-?significant|"
    r"bonferroni|holm|tukey|benjamini|FDR|false discovery|post[- ]hoc|emmeans|multcomp|"
    r"covariates?|predictors?|fixed effects?|random effects?|sample size|n ?= ?\d+)\b")


def to_text(path):
    path = Path(path)
    txt = path.with_suffix(".txt")
    if not txt.exists():
        if path.suffix == ".pdf":
            subprocess.run(["pdftotext", str(path), str(txt)], check=True)
        else:
            x = path.read_text()
            body = x[x.find("<body"):x.rfind("</body>")]
            body = re.sub(r"<title>", "\n## ", body)
            t = html.unescape(re.sub(r"<[^>]+>", " ", body))
            txt.write_text(re.sub(r"[ \t]+", " ", t))
    return txt.read_text(errors="replace")


def main():
    path = sys.argv[1]
    max_chars = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
    t = to_text(path)
    ref = re.search(r"\n\s*(References|REFERENCES|Literature Cited|LITERATURE CITED)\s*\n", t)
    if ref:
        t = t[:ref.start()]
    flat = re.sub(r"-\n(?=[a-z])", "", t)
    flat = re.sub(r"\s+", " ", flat)

    m = None
    for pat in (r"(?i)statistical analys[ie]s\W", r"(?i)data analys[ie]s\W",
                r"(?i)statistical methods\W", r"(?i)\bstatistics\W(?!\s*,)"):
        m = re.search(pat, flat)
        if m:
            break
    if m:
        print("=== STATS SECTION ===")
        print(flat[m.start():m.start() + max_chars])
    print("\n=== KEYWORD SENTENCES ===")
    seen = 0
    for s in re.split(r"(?<=[.!?])\s+(?=[A-Z])", flat):
        if KEYWORDS.search(s) and len(s) < 600:
            print("-", s)
            seen += len(s)
            if seen > max_chars:
                break
    print(f"\n[{len(flat.split())} words]")


if __name__ == "__main__":
    main()
