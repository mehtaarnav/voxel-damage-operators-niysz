"""Environment check: report presence and version of every package the study needs."""
import importlib

MODS = [
    "numpy", "scipy", "skimage", "networkx", "tifffile",
    "skan", "matplotlib", "pandas", "openpyxl", "requests",
]

for m in MODS:
    try:
        mod = importlib.import_module(m)
        print("{:12s} OK       {}".format(m, getattr(mod, "__version__", "?")))
    except Exception as e:
        print("{:12s} MISSING  ({}: {})".format(m, type(e).__name__, e))
