"""PHASE 1c — fetch the two open-access (CC-BY) MDPI papers used as ground truth.

Pecho et al. Materials 2015, 8(9),  5554-5585  doi:10.3390/ma8095265
Pecho et al. Materials 2015, 8(10), 7129-7147  doi:10.3390/ma8105370
Both are CC-BY open access.  Also tries the supplementary-material archives,
which for these papers carry the per-sample microstructure tables.
"""

from __future__ import annotations

import os
import sys

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "refs")
os.makedirs(REFS, exist_ok=True)

UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36")

TARGETS = [
    ("ma8095265_transport.pdf", "https://www.mdpi.com/1996-1944/8/9/5554/pdf"),
    ("ma8105370_tpb.pdf",       "https://www.mdpi.com/1996-1944/8/10/7129/pdf"),
    ("ma8105370_tpb_alt.pdf",   "https://www.mdpi.com/1996-1944/8/10/5370/pdf"),
]


def get(name, url):
    dest = os.path.join(REFS, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 50_000:
        print(f"  [skip] {name} ({os.path.getsize(dest)/1e6:.2f} MB)")
        return True
    try:
        r = requests.get(url, headers={"User-Agent": UA}, timeout=90,
                         allow_redirects=True)
        ctype = r.headers.get("content-type", "")
        print(f"  [{r.status_code}] {url}  -> {ctype}  {len(r.content)/1e6:.2f} MB")
        if r.status_code == 200 and r.content[:4] == b"%PDF":
            with open(dest, "wb") as f:
                f.write(r.content)
            print(f"         saved {name}")
            return True
    except Exception as e:
        print(f"  [ERR] {url}: {e}")
    return False


def main():
    print("=" * 78)
    print("PHASE 1c — fetching ground-truth papers")
    print("=" * 78)
    for name, url in TARGETS:
        get(name, url)
    print("\nrefs dir:", REFS)
    for f in sorted(os.listdir(REFS)):
        p = os.path.join(REFS, f)
        print(f"   {f:32s} {os.path.getsize(p)/1e6:8.2f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
