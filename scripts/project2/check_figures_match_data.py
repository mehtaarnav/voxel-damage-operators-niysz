"""Check that committed figures still show what their source data says.

The obvious check -- regenerate and compare bytes -- does not work across
platforms. matplotlib embeds glyph outlines produced by FreeType, and different
FreeType versions produce different outlines for identical content, so a figure
built on Windows never matches one built on Linux byte-for-byte. matplotlib's
own test suite pins FreeType for this reason; pip cannot.

Comparing rendered pixels with a tolerance tests the thing that actually
matters. A figure regenerated from changed data moves marks, which shows up as
a large difference. A figure regenerated on a different FreeType shifts glyph
edges by a fraction of a pixel, which does not.

Usage: python check_figures_match_data.py [--update]
  default    regenerate into a temporary directory and compare against the
             committed figures, failing if any differs beyond tolerance
  --update   overwrite the committed figures (use after a deliberate change)
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIGS = os.path.join(ROOT, "out", "writeup", "figs")
SCRIPTS = ["make_figures.py", "make_figures_o7.py", "make_figure_threshold.py"]

# mean absolute pixel difference, 0-255 scale. Glyph-edge differences from a
# different FreeType land around 0.1-0.5; a moved data mark lands far above.
TOL = 2.0


def raster(pdf, out_prefix, dpi=80):
    subprocess.run(["pdftoppm", "-r", str(dpi), "-png", "-singlefile", pdf,
                    out_prefix], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return out_prefix + ".png"


def load(png):
    from PIL import Image
    return np.asarray(Image.open(png).convert("L"), dtype=float)


def main():
    if "--update" in sys.argv:
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            [ROOT, os.path.join(ROOT, "scripts", "project2")]),
            MPLBACKEND="Agg")
        for s in SCRIPTS:
            subprocess.run([sys.executable, os.path.join(ROOT, "scripts",
                                                         "project2", s)],
                           check=True, env=env)
        print("figures regenerated in place")
        return 0

    tmp = tempfile.mkdtemp()
    backup = os.path.join(tmp, "committed")
    shutil.copytree(FIGS, backup)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(
        [ROOT, os.path.join(ROOT, "scripts", "project2")]), MPLBACKEND="Agg")
    try:
        for s in SCRIPTS:
            subprocess.run([sys.executable, os.path.join(ROOT, "scripts",
                                                         "project2", s)],
                           check=True, env=env,
                           stdout=subprocess.DEVNULL)
        bad = []
        names = sorted(f for f in os.listdir(backup) if f.endswith(".pdf"))
        print(f"{'figure':<32}{'mean |diff|':>12}")
        print("-" * 44)
        for n in names:
            a = raster(os.path.join(backup, n), os.path.join(tmp, "a"))
            b = raster(os.path.join(FIGS, n), os.path.join(tmp, "b"))
            ia, ib = load(a), load(b)
            if ia.shape != ib.shape:
                d = float("inf")
            else:
                d = float(np.abs(ia - ib).mean())
            flag = "" if d <= TOL else "   <- DIFFERS"
            if d > TOL:
                bad.append(n)
            print(f"{n:<32}{d:>12.3f}{flag}")
        print()
        if bad:
            print("These figures no longer match their source data:")
            for n in bad:
                print("   ", n)
            print("\nRegenerate with --update and commit, or fix the script.")
            return 1
        print(f"all {len(names)} figures match within tolerance {TOL}")
        return 0
    finally:
        # restore the committed versions so a check never dirties the tree
        for f in os.listdir(backup):
            shutil.copy2(os.path.join(backup, f), os.path.join(FIGS, f))
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
