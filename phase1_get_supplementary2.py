"""PHASE 1e(ii) — second attempt at the supplementary tables, via Europe PMC.

The NCBI oa_package tar.gz 404s for both PMCIDs.  Europe PMC exposes a
supplementaryFiles endpoint that returns a zip of the article's supplementary
material, which is where Tables S1/S2 (the per-sample numeric tables) live.
"""

from __future__ import annotations

import io
import os
import sys
import zipfile

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "refs")
SUPP = os.path.join(REFS, "supplementary")
os.makedirs(SUPP, exist_ok=True)

UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36")

PMCIDS = {"PMC5455394": "ma8105370_TPB", "PMC5512617": "ma8095265_transport"}

EPMC_SUPP = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
EPMC_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"


def try_supp(pmcid, tag):
    url = EPMC_SUPP.format(pmcid=pmcid)
    try:
        r = requests.get(url, headers={"User-Agent": UA}, timeout=180)
    except Exception as e:
        print(f"   [ERR] {e}")
        return False
    print(f"   supplementaryFiles [{r.status_code}]  {len(r.content)/1e3:.1f} kB  "
          f"{r.headers.get('content-type')}")
    if r.status_code != 200 or len(r.content) < 1000:
        return False
    outdir = os.path.join(SUPP, tag)
    os.makedirs(outdir, exist_ok=True)
    try:
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            for n in z.namelist():
                print(f"      member: {n}  ({z.getinfo(n).file_size/1e3:.1f} kB)")
            z.extractall(outdir)
        return True
    except zipfile.BadZipFile:
        dest = os.path.join(outdir, f"{pmcid}_supp.bin")
        with open(dest, "wb") as f:
            f.write(r.content)
        print(f"      not a zip; saved raw to {dest}")
        return False


def try_fulltext(pmcid, tag):
    url = EPMC_XML.format(pmcid=pmcid)
    r = requests.get(url, headers={"User-Agent": UA}, timeout=120)
    print(f"   fullTextXML [{r.status_code}]  {len(r.content)/1e3:.1f} kB")
    if r.status_code == 200 and len(r.content) > 10_000:
        dest = os.path.join(REFS, f"{tag}_{pmcid}.epmc.xml")
        with open(dest, "wb") as f:
            f.write(r.content)
        print(f"      saved {os.path.basename(dest)}")


def main():
    print("=" * 78)
    print("PHASE 1e(ii) — Europe PMC supplementary + full text")
    print("=" * 78)
    for pmcid, tag in PMCIDS.items():
        print(f"\n{pmcid} ({tag})")
        try_supp(pmcid, tag)
        try_fulltext(pmcid, tag)

    print("\nSupplementary tree:")
    n = 0
    for dp, _, fns in os.walk(SUPP):
        for fn in fns:
            p = os.path.join(dp, fn)
            print(f"   {os.path.relpath(p, SUPP):60s} {os.path.getsize(p)/1e3:9.1f} kB")
            n += 1
    if n == 0:
        print("   (empty)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
