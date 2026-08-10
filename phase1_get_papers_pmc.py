"""PHASE 1c(ii) — pull RAW full text (incl. tables) of the two ground-truth papers.

MDPI is behind Cloudflare, so we use NCBI PMC, where both papers are deposited
open-access:

  PMC5455394  Pecho et al., Materials 2015, 8(10), 7129-7147
              "3D Microstructure Effects in Ni-YSZ Anodes: Influence of TPB
               Lengths on the Electrochemical Performance"   doi:10.3390/ma8105370
  PMC5512617  Pecho et al., Materials 2015, 8(9), 5554-5585
              "3D Microstructure Effects in Ni-YSZ Anodes: Prediction of Effective
               Transport Properties and Optimization of Redox Stability"
                                                             doi:10.3390/ma8095265

We deliberately fetch the machine-readable BioC XML rather than scraping the
rendered HTML, and rather than trusting an LLM summary of the page: table cells
must be read verbatim because they are the ground truth for the Phase-2 and
Phase-4 gates.
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

PAPERS = {
    "PMC5455394": "ma8105370_TPB",
    "PMC5512617": "ma8095265_transport",
}

ENDPOINTS = [
    ("bioc.xml",
     "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmcid}/unicode"),
    ("efetch.xml",
     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid}&retmode=xml"),
]


def fetch(pmcid, tag):
    got_any = False
    for suffix, tmpl in ENDPOINTS:
        url = tmpl.format(pmcid=pmcid)
        dest = os.path.join(REFS, f"{tag}_{pmcid}.{suffix}")
        if os.path.exists(dest) and os.path.getsize(dest) > 20_000:
            print(f"  [skip] {os.path.basename(dest)} "
                  f"({os.path.getsize(dest)/1e3:.0f} kB)")
            got_any = True
            continue
        try:
            r = requests.get(url, headers={"User-Agent": UA}, timeout=120)
            print(f"  [{r.status_code}] {url[:80]}...  {len(r.content)/1e3:.0f} kB")
            if r.status_code == 200 and len(r.content) > 5_000:
                with open(dest, "wb") as f:
                    f.write(r.content)
                print(f"         saved {os.path.basename(dest)}")
                got_any = True
        except Exception as e:
            print(f"  [ERR] {url}: {e}")
    return got_any


def main():
    print("=" * 78)
    print("PHASE 1c(ii) — raw full text from PMC")
    print("=" * 78)
    for pmcid, tag in PAPERS.items():
        print(f"\n{pmcid}  ({tag})")
        fetch(pmcid, tag)
    print("\nrefs dir contents:")
    for f in sorted(os.listdir(REFS)):
        p = os.path.join(REFS, f)
        print(f"   {f:44s} {os.path.getsize(p)/1e3:9.1f} kB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
