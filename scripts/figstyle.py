"""Shared APS-style matplotlib configuration."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito (colorblind-safe), fixed assignment order
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7",
     "#56B4E9", "#333333"]

plt.rcParams.update({
    "font.size": 7.5,
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "axes.linewidth": 0.6,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
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

DBL = 7.05   # double-column width (in)
SGL = 3.40   # single-column width


def panel_label(ax, s, dx=-0.14, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left")
