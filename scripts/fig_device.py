"""Figure 1: three-dimensional device schematic of the graphene
Josephson junction thermal detector.

Every labeled fact is a published parameter of the Ta/Ti/Au recipe of
Jung et al. (Table I): Ta(10 nm)/Ti(60 nm)/Au(5 nm) contacts on
monolayer graphene, channel L = 0.2 um, W = 5.3 um, on 280 nm SiO2
with a doped-Si back gate; Delta* = 1.764 kB Tc* = 87 ueV for
Tc* = 0.57 K. The readout embedding (6 GHz resonator, 2 nH external
inductance) and the electron-phonon cooling law Sigma A (Te^3 - Tp^3)
are the ones used throughout the manuscript. Geometry is drawn not to
scale.

Only exposed surfaces are drawn (the SiO2 top is a frame around the
graphene footprint, the graphene top only where the contacts do not
cover it), so the depth sort never has to resolve hidden coplanar
faces. The 2D annotation overlay projects anchor points with the
finalized 3D transform (fig.canvas.draw() precedes the projection).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import proj3d
from figstyle import C, SGL
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

FS = 5.8
LC = "0.3"

GX = (0.5, 9.5); GY = (0.8, 6.2)          # graphene sheet footprint
XL = (0.7, 4.15); XR = (5.85, 9.3)        # contact x-extents
YC = (1.0, 6.0)                           # contact y-extent
ZG = 0.07                                 # graphene (exaggerated)
ZTA, ZTI, ZAU = 0.22, 0.80, 0.92          # cumulative layer tops

verts, cols = [], []


def shade(c, f):
    r, g, b = mcolors.to_rgb(c)
    return (min(r * f, 1), min(g * f, 1), min(b * f, 1))


def quad(p0, p1, p2c, p3, c):
    verts.append([p0, p1, p2c, p3]); cols.append(c)


def top(xr, yr, z, c):
    quad((xr[0], yr[0], z), (xr[1], yr[0], z),
         (xr[1], yr[1], z), (xr[0], yr[1], z), shade(c, 1.0))


def front(xr, y, zr, c):
    quad((xr[0], y, zr[0]), (xr[1], y, zr[0]),
         (xr[1], y, zr[1]), (xr[0], y, zr[1]), shade(c, 0.80))


def right(x, yr, zr, c):
    quad((x, yr[0], zr[0]), (x, yr[1], zr[0]),
         (x, yr[1], zr[1]), (x, yr[0], zr[1]), shade(c, 0.62))


def left(x, yr, zr, c):
    quad((x, yr[0], zr[0]), (x, yr[1], zr[0]),
         (x, yr[1], zr[1]), (x, yr[0], zr[1]), shade(c, 0.72))


SI, OX, GR = "#9aa1ab", "#d3dfeb", "#3a3a3a"
TA, TI, AU = "#8fa0b8", "#bcc5cc", "#e6bf55"

# doped Si substrate (its top is fully covered by the SiO2)
front((0, 10), 0.0, (-1.45, -0.72), SI)
right(10.0, (0, 7), (-1.45, -0.72), SI)
left(0.0, (0, 7), (-1.45, -0.72), SI)
# SiO2: sides plus a top frame around the graphene footprint
front((0, 10), 0.0, (-0.72, 0.0), OX)
right(10.0, (0, 7), (-0.72, 0.0), OX)
left(0.0, (0, 7), (-0.72, 0.0), OX)
top((0, 10), (0, GY[0]), 0.0, OX)
top((0, 10), (GY[1], 7), 0.0, OX)
top((0, GX[0]), GY, 0.0, OX)
top((GX[1], 10), GY, 0.0, OX)
# graphene: outer faces plus the exposed top pieces
front(GX, GY[0], (0.0, ZG), GR)
right(GX[1], GY, (0.0, ZG), GR)
left(GX[0], GY, (0.0, ZG), GR)
top((XL[1], XR[0]), (GY[0], GY[1]), ZG, GR)        # junction channel
top((GX[0], XL[0]), GY, ZG, GR)
top((XR[1], GX[1]), GY, ZG, GR)
top((XL[0], XR[1]), (GY[0], YC[0]), ZG, GR)
top((XL[0], XR[1]), (YC[1], GY[1]), ZG, GR)
# contacts: Ta / Ti / Au stacks
for xr in (XL, XR):
    for z0, z1, c in ((ZG, ZTA, TA), (ZTA, ZTI, TI), (ZTI, ZAU, AU)):
        front(xr, YC[0], (z0, z1), c)
        right(xr[1], YC, (z0, z1), c)
        left(xr[0], YC, (z0, z1), c)
    top(xr, YC, ZAU, AU)

fig = plt.figure(figsize=(SGL, 2.95))
under = fig.add_axes([0, 0, 1, 1]); under.axis("off")
under.set_xlim(0, 1); under.set_ylim(0, 1)
ax = fig.add_axes([-0.18, 0.06, 0.90, 0.86], projection="3d")
ax.set_axis_off()
ax.view_init(elev=18, azim=-70)
ax.set_proj_type("ortho")
ax.set_box_aspect((10, 7, 3.6), zoom=1.00)
pc = Poly3DCollection(verts, facecolors=cols, edgecolors="0.35",
                      linewidths=0.3, zsort="average")
ax.add_collection3d(pc)

# honeycomb lattice hint on the exposed channel (view-culled)
a = 0.16
gx0, gx1 = XL[1] + 0.10, XR[0] - 0.10
gy0, gy1 = YC[0] + 0.10, YC[1] - 0.10


def visible(p):
    if not (gx0 <= p[0] <= gx1 and gy0 <= p[1] <= gy1):
        return False
    return not (p[0] > 4.9 and p[1] > 0.7 + 2.747 * (XR[0] - p[0]))


for iy in range(int((gy1 - gy0) / (a * np.sqrt(3))) + 2):
    for ix in range(int((gx1 - gx0) / (1.5 * a)) + 2):
        cx = gx0 + 1.5 * a * ix
        cy = gy0 + a * np.sqrt(3) * (iy + 0.5 * (ix % 2))
        pts = [(cx + a * np.cos(np.pi / 3 * k),
                cy + a * np.sin(np.pi / 3 * k)) for k in range(7)]
        for p0, p1 in zip(pts[:-1], pts[1:]):
            if visible(p0) and visible(p1):
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [ZG + 0.02] * 2,
                        color="0.66", lw=0.3, zorder=50)

ax.set_xlim(0, 10); ax.set_ylim(0, 7); ax.set_zlim(-1.5, 2.6)
fig.canvas.draw()          # finalize the 3D transform before projecting

# soft drop shadow under the substrate (underlay, beneath the 3D axes)
def pshadow(x, y, z):
    tx, ty, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
    return tuple(fig.transFigure.inverted().transform(
        ax.transData.transform((tx, ty))))


base = [pshadow(*pt) for pt in ((0, 0, -1.45), (10, 0, -1.45),
                                (10, 7, -1.45), (0, 7, -1.45))]
for off, al in (((0.006, -0.010), 0.10), ((0.013, -0.020), 0.06)):
    under.fill([b[0] + off[0] for b in base],
               [b[1] + off[1] for b in base], color="k", alpha=al,
               lw=0)

# ---- 2D annotation overlay -------------------------------------------
ov = fig.add_axes([0, 0, 1, 1]); ov.axis("off")
ov.set_xlim(0, 1); ov.set_ylim(0, 1)


def p2(x, y, z):
    tx, ty, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
    return tuple(fig.transFigure.inverted().transform(
        ax.transData.transform((tx, ty))))


def facetext(pt3, s, fs=FS, color="0.12"):
    x, y = p2(*pt3)
    ov.text(x, y, s, fontsize=fs, color=color, ha="center",
            va="center")


# phases and gate voltage, written directly on the layers ---------------
facetext((2.9, 3.2, ZAU), r"$-\varphi/2$", fs=6.2)
facetext((8.1, 2.4, ZAU), r"$+\varphi/2$", fs=6.2)
facetext((5.0, 0.0, -1.13), r"$V_{\rm BG}=30$ V", fs=6.0,
         color="0.97")

# process arrows on the device (text lives in the key) ------------------
q0, q1 = p2(3.55, 3.6, 0.34), p2(6.45, 3.6, 0.34)
ov.annotate("", xy=q1, xytext=q0,
            arrowprops=dict(arrowstyle="-|>", lw=1.1, color=C[1],
                            mutation_scale=9))
e0, e1 = p2(5.45, 5.2, 0.12), p2(7.25, 5.2, 0.72)
ov.annotate("", xy=e1, xytext=e0,
            arrowprops=dict(arrowstyle="<|-|>", lw=0.9, color=C[0],
                            mutation_scale=7))
c0, c1 = p2(6.55, 0.9, 0.06), p2(6.55, 0.9, -1.30)
ov.annotate("", xy=c1, xytext=c0,
            arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C[2],
                            mutation_scale=7))


def squiggle(start, end, amp=0.011, cyc=6.5, col=C[4], lw=0.9):
    tt = np.linspace(0, 1, 160)
    wob = amp * np.sin(2 * np.pi * cyc * tt) * (1 - 0.65 * tt)
    dxp, dyp = end[0] - start[0], end[1] - start[1]
    nn = np.hypot(dxp, dyp)
    ux, uy = dxp / nn, dyp / nn
    px = start[0] + dxp * tt - uy * wob
    py = start[1] + dyp * tt + ux * wob
    ov.plot(px[:-8], py[:-8], color=col, lw=lw)
    ov.annotate("", xy=end, xytext=(px[-9], py[-9]),
                arrowprops=dict(arrowstyle="-|>", lw=lw, color=col,
                                mutation_scale=7))


squiggle((0.415, 0.825), p2(4.50, 5.05, 0.11), cyc=5)

# readout resonator (top left, wired to the left contact) ---------------
r0 = p2(2.4, 5.6, ZAU)
rx, ry = 0.170, 0.930
ov.plot([r0[0], rx - 0.055], [r0[1], ry], color="0.15", lw=0.7)
th = np.linspace(np.pi, 0, 40)
for k in range(4):
    ov.plot(rx - 0.049 + 0.026 * k + 0.013 * np.cos(th),
            ry + 0.016 * np.sin(th), color="0.15", lw=0.7)
ov.plot([rx + 0.055, rx + 0.072], [ry, ry], color="0.15", lw=0.7)
for xc in (rx + 0.072, rx + 0.084):
    ov.plot([xc, xc], [ry - 0.016, ry + 0.016], color="0.15", lw=0.7)
ov.plot([rx + 0.084, rx + 0.100], [ry, ry], color="0.15", lw=0.7)
ov.text(rx + 0.115, ry,
        "readout resonator,\n$\\nu_r=6$ GHz, $L_r=2$ nH",
        fontsize=FS, ha="left", va="center")

ov.text(0.990, 0.212, "not to scale", fontsize=4.8,
        color="0.45", ha="right")

# ---- key (legend-style labeling, no leader lines) ---------------------
from matplotlib.patches import FancyBboxPatch, Rectangle

LX, LT = 0.678, 0.712            # swatch x, text x
FL = 5.6                         # legend font size
box = FancyBboxPatch((0.663, 0.225), 0.330, 0.730,
                     boxstyle="round,pad=0.008,rounding_size=0.012",
                     fc="white", ec="0.80", lw=0.6,
                     transform=ov.transAxes, zorder=4)
ov.add_patch(box)


def swatch(y, c, label):
    ov.add_patch(Rectangle((LX, y - 0.011), 0.022, 0.022, fc=c,
                           ec="0.35", lw=0.4, transform=ov.transAxes,
                           zorder=5))
    ov.text(LT, y, label, fontsize=FL, va="center", zorder=5)


ov.text(LX, 0.925, "S contact,  $\\Delta^*=87\\ \\mu$eV:",
        fontsize=FL, va="center", color="0.15", zorder=5)
swatch(0.878, AU, "Au (5 nm)")
swatch(0.831, TI, "Ti (60 nm)")
swatch(0.784, TA, "Ta (10 nm, adhesion)")
swatch(0.722, GR, "monolayer graphene")
swatch(0.675, OX, r"SiO$_2$ (280 nm)")
swatch(0.628, SI, "doped Si back gate")
ov.text(LX, 0.570, r"$L=0.2\ \mu$m,   $W=5.3\ \mu$m",
        fontsize=FL, va="center", color="0.15", zorder=5)


def glyphrow(y, kind, col, label, dy=0.0):
    xa, xb = LX - 0.008, LX + 0.046
    if kind == "arrow":
        ov.annotate("", xy=(xb, y), xytext=(xa, y), zorder=5,
                    arrowprops=dict(arrowstyle="-|>", lw=0.9,
                                    color=col, mutation_scale=7))
    elif kind == "double":
        ov.annotate("", xy=(xb, y), xytext=(xa, y), zorder=5,
                    arrowprops=dict(arrowstyle="<|-|>", lw=0.9,
                                    color=col, mutation_scale=6))
    elif kind == "wavy":
        tt = np.linspace(0, 1, 80)
        ov.plot(xa + (xb - xa) * tt,
                y + 0.006 * np.sin(2 * np.pi * 3 * tt), color=col,
                lw=0.9, zorder=5)
    ov.text(LT + 0.020, y + dy, label, fontsize=FL, va="center",
            zorder=5)


glyphrow(0.500, "arrow", C[1], r"supercurrent $I_s(\varphi,T_e)$")
glyphrow(0.436, "double", C[0],
         "quasiparticle exchange,\n$\\tau_{\\rm A}$", dy=-0.012)
glyphrow(0.352, "wavy", C[4], "microwave photon, $h\\nu$")
glyphrow(0.288, "arrow", C[2],
         "e-ph cooling,\n$\\Sigma A\\,(T_e^3-T_p^3)$", dy=-0.012)

# ---- Andreev-doublet energy inset (bottom left) -----------------------
ibox = FancyBboxPatch((0.030, 0.045), 0.295, 0.225,
                      boxstyle="round,pad=0.008,rounding_size=0.012",
                      fc="white", ec="0.80", lw=0.6,
                      transform=ov.transAxes, zorder=4)
ov.add_patch(ibox)
ov.text(0.190, 0.252, "Andreev doublet (per channel)", fontsize=5.2,
        ha="center", va="center", color="0.15", zorder=5)
bx0, bx1 = 0.078, 0.215
Etop, Ebot, Ec = 0.205, 0.145, 0.175
ov.add_patch(Rectangle((bx0, Etop), bx1 - bx0, 0.030, fc="0.88",
                       ec="none", transform=ov.transAxes, zorder=5))
ov.add_patch(Rectangle((bx0, Ebot - 0.030), bx1 - bx0, 0.030,
                       fc="0.88", ec="none", transform=ov.transAxes,
                       zorder=5))
ov.text(0.5 * (bx0 + bx1), Etop + 0.015, "continuum", fontsize=4.8,
        ha="center", va="center", color="0.40", zorder=6)
ov.text(0.5 * (bx0 + bx1), Ebot - 0.015, "continuum", fontsize=4.8,
        ha="center", va="center", color="0.40", zorder=6)
for yy in (Etop, Ebot):
    ov.plot([bx0, bx1], [yy, yy], color="0.35", lw=0.6, ls=(0, (3, 2)),
            zorder=6)
for sgn in (+1, -1):
    ov.plot([bx0 + 0.010, bx1 - 0.010],
            [Ec + sgn * 0.019] * 2, color=C[0], lw=1.4, zorder=6)
ov.text(0.223, Etop, r"$+\Delta^*$", fontsize=5.2, va="center",
        zorder=5)
ov.text(0.223, Ec, r"$\pm E_m(\varphi)$", fontsize=5.2, va="center",
        color=C[0], zorder=5)
ov.text(0.223, Ebot, r"$-\Delta^*$", fontsize=5.2, va="center",
        zorder=5)
ov.annotate("", xy=(0.066, 0.230), xytext=(0.066, 0.120), zorder=5,
            arrowprops=dict(arrowstyle="-|>", lw=0.6, color="0.30",
                            mutation_scale=5))
ov.text(0.057, 0.220, "$E$", fontsize=5.2, ha="right", va="center",
        zorder=5)
ov.text(0.178, 0.072,
        r"$E_m=\pm\Delta^*\sqrt{1-\tau\sin^2(\varphi/2)}$",
        fontsize=5.3, ha="center", va="center", zorder=5)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig_device.pdf"))
print("fig_device done")
