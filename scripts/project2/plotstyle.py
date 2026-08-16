"""Shared plotting style for the Project 2 manuscript figures.

Serif text to match the elsarticle body font, hairline spines, outward minor
ticks, no grid. The three-anode palette is a monotone blue ramp rather than a
categorical scheme, because coarseness is an ordinal variable and the ramp
survives greyscale printing. Synthetic-platform data is always grey, so that
"real" versus "synthetic" reads before the legend does.
"""
import matplotlib as mpl

ANODES = ["fine", "medium", "coarse"]
COLOR = {"fine": "#08306B", "medium": "#2171B5", "coarse": "#6BAED6"}
MARKER = {"fine": "o", "medium": "s", "coarse": "^"}
SYNTH = "#8C8C8C"
ACCENT = "#B2182B"          # violations, forbidden regions, measured targets
COL_W = 3.35                # single column, inches (Elsevier 5.5 pica ~ 88 mm)
FULL_W = 6.95               # double column, inches (~176 mm)

RC = {
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "axes.grid": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlelocation": "left",
    "axes.titlepad": 4.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.45,
    "ytick.minor.width": 0.45,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.6,
    "ytick.minor.size": 1.6,
    "lines.linewidth": 1.1,
    "lines.markersize": 3.6,
    "legend.frameon": False,
    "legend.handlelength": 1.6,
    "legend.labelspacing": 0.3,
    "figure.dpi": 400,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
}


def use():
    mpl.rcParams.update(RC)


def panel(ax, letter, dx=-0.155, dy=1.06):
    """Panel letter, placed on the title's own baseline.

    The letter used to be a free-floating text at a fixed offset in axes
    coordinates while the title was positioned by matplotlib's title padding.
    Those are two independent placements, so whether they appeared level
    depended on the figure's height, and across a paper with panels of
    different heights they visibly did not.

    Writing the letter into the title guarantees a shared baseline. Bold comes
    from mathtext so the rest of the title keeps its normal weight.
    """
    # NB: rcParams sets titlelocation to "left", so set_title() writes the LEFT
    # title slot. A bare get_title() reads the CENTRE slot and would return an
    # empty string, silently discarding the panel's real title.
    existing = ax.get_title(loc="left")
    ax.set_title(rf"$\bf{{({letter})}}$   {existing}".rstrip(), loc="left")


def tidy(ax, minor_x=True, minor_y=True):
    from matplotlib.ticker import AutoMinorLocator
    if minor_x and ax.get_xscale() == "linear":
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    if minor_y and ax.get_yscale() == "linear":
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
