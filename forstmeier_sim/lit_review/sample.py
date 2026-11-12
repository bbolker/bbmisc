"""Draw the coding sample from screen.csv: 10 empirical papers per journal at
random (all of them when a journal has fewer than 10). Writes sample.csv with
every empirical paper in seeded random order per journal; `sampled` marks the
first 10. If a sampled paper turns out not to be empirical on reading, mark it
excluded and take the next paper in `rank` order for that journal.
"""

import csv
import random
from pathlib import Path

HERE = Path(__file__).parent
SEED = 20261112
N_PER_JOURNAL = 10

rows = [r for r in csv.DictReader(open(HERE / "screen.csv")) if r["empirical"] == "y"]
rng = random.Random(SEED)

out = []
for journal in sorted({r["journal"] for r in rows}):
    sub = sorted((r for r in rows if r["journal"] == journal), key=lambda r: r["doi"])
    rng.shuffle(sub)
    for rank, r in enumerate(sub, start=1):
        out.append({"journal": journal, "rank": rank, "sampled": "y" if rank <= N_PER_JOURNAL else "n",
                    "doi": r["doi"], "title": r["title"]})

with open(HERE / "sample.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)

for journal in sorted({o["journal"] for o in out}):
    n = sum(o["journal"] == journal and o["sampled"] == "y" for o in out)
    print(f"{journal:10s} {n} sampled")
