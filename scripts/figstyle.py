"""Shared APS-style matplotlib configuration."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito (colorblind-safe), fixed assignment order
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7",
     "#56B4E9", "#333333"]

plt.rcParams.update({
    "font.size": 8,
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "axes.linewidth": 0.6,
    "axes.labelsize": 8.5,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "legend.frameon": False,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "lines.linewidth": 1.1,
    "savefig.dpi": 600,
    "figure.dpi": 130,
})

# Figures are drawn at their printed size: the text width of the
# manuscript (elsarticle preprint, 390 pt = 5.40 in), so that the font
# sizes above are the sizes on the page.
DBL = 5.40   # full text width (in)
SGL = 3.40   # narrow figure
FS_IN = 7.0  # minimum size of in-panel annotations (pt)


def panel_label(ax, s, dx=-0.16, dy=1.05):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="top", ha="left")
