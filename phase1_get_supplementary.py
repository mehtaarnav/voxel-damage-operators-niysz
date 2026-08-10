"""PHASE 1e — retrieve the supplementary material (the numeric tables).

Neither paper's MAIN TEXT carries per-sample numeric tables: the transport paper
(ma8095265) has exactly one table and it is qualitative (+/- symbols), and the
TPB paper (ma8105370) has no <table-wrap> at all.  The numbers we need for the
ground-truth table live in:
    Table S1, S2  (ma8105370 supplementary, id materials-08-05370-s001)
    supplementary (ma8095265, id materials-08-05265-s001)

The PMC Open Access service exposes the complete article package (including
supplementary files) as a tar.gz.  We use that rather than MDPI, which is behind
Cloudflare.
"""

from __future__ import annotations

import os
import sys
import tarfile
import xml.etree.ElementTree as ET

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "refs")
SUPP = os.path.join(REFS, "supplementary")
os.makedirs(SUPP, exist_ok=True)

UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36")

PMCIDS = {"PMC5455394": "ma8105370_TPB", "PMC5512617": "ma8095265_transport"}

OA_API = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}"


def oa_links(pmcid):
    r = requests.get(OA_API.format(pmcid=pmcid), headers={"User-Agent": UA},
                     timeout=60)
    print(f"  oa.fcgi [{r.status_code}] for {pmcid}")
    if r.status_code != 200:
        return []
    try:
        root = ET.fromstring(r.content)
    except ET.ParseError as e:
        print("   parse error:", e)
        print(r.text[:500])
        return []
    err = root.find(".//error")
    if err is not None:
        print("   OA service error:", err.attrib, err.text)
        return []
    out = []
    for link in root.iter("link"):
        out.append((link.attrib.get("format"), link.attrib.get("href")))
        print(f"   link format={link.attrib.get('format')}  "
              f"href={link.attrib.get('href')}")
    return out


def fetch_package(pmcid, tag, links):
    tgz = [h for f, h in links if f == "tgz"]
    if not tgz:
        print("   no tgz package offered")
        return False
    url = tgz[0].replace("ftp://ftp.ncbi.nlm.nih.gov",
                         "https://ftp.ncbi.nlm.nih.gov")
    dest = os.path.join(REFS, f"{tag}_{pmcid}.tar.gz")
    if not (os.path.exists(dest) and os.path.getsize(dest) > 10_000):
        r = requests.get(url, headers={"User-Agent": UA}, timeout=300, stream=True)
        print(f"   package [{r.status_code}] {url}")
        if r.status_code != 200:
            return False
        with open(dest, "wb") as f:
            for c in r.iter_content(1 << 20):
                f.write(c)
    print(f"   package {os.path.getsize(dest)/1e6:.2f} MB")

    outdir = os.path.join(SUPP, tag)
    os.makedirs(outdir, exist_ok=True)
    with tarfile.open(dest, "r:gz") as t:
        members = t.getmembers()
        print(f"   package contains {len(members)} members:")
        for m in members:
            print(f"      {m.size/1e3:9.1f} kB  {m.name}")
        t.extractall(outdir, filter="data")
    return True


def main():
    print("=" * 78)
    print("PHASE 1e — supplementary material retrieval")
    print("=" * 78)
    for pmcid, tag in PMCIDS.items():
        print(f"\n{pmcid} ({tag})")
        links = oa_links(pmcid)
        if links:
            fetch_package(pmcid, tag, links)
    print("\nSupplementary tree:")
    for dp, _, fns in os.walk(SUPP):
        for fn in fns:
            p = os.path.join(dp, fn)
            print(f"   {os.path.relpath(p, SUPP):60s} {os.path.getsize(p)/1e3:9.1f} kB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
