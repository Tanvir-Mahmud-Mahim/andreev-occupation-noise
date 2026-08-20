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
ax = fig.add_axes([-0.13, -0.02, 1.06, 0.86], projection="3d")
ax.set_axis_off()
ax.view_init(elev=18, azim=-70)
ax.set_proj_type("ortho")
ax.set_box_aspect((10, 7, 3.6), zoom=1.12)
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

# ---- 2D annotation overlay -------------------------------------------
ov = fig.add_axes([0, 0, 1, 1]); ov.axis("off")
ov.set_xlim(0, 1); ov.set_ylim(0, 1)


def p2(x, y, z):
    tx, ty, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
    return tuple(fig.transFigure.inverted().transform(
        ax.transData.transform((tx, ty))))


def callout(text, pt3, txt2, ha="left", fs=FS, color="0.1"):
    ov.annotate(text, xy=p2(*pt3), xytext=txt2, ha=ha, va="center",
                fontsize=fs, color=color,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=LC,
                                shrinkA=1.5, shrinkB=0.5))


def facetext(pt3, s, fs=FS, color="0.12", dx=0.0, dy=0.0, rot=0.0):
    x, y = p2(*pt3)
    ov.text(x + dx, y + dy, s, fontsize=fs, color=color, ha="center",
            va="center", rotation=rot, rotation_mode="anchor")


# labels written directly on the layers ---------------------------------
facetext((2.4, 3.9, ZAU), "S  (Ta/Ti/Au)\n$\\Delta^*=87\\ \\mu$eV",
         fs=6.0)
facetext((3.6, 1.3, ZAU), r"$-\varphi/2$", fs=6.2)
facetext((8.0, 2.3, ZAU), r"$+\varphi/2$", fs=6.2)
facetext((1.8, 0.0, -0.36), r"SiO$_2$ (280 nm)", fs=6.0)
facetext((4.0, 0.0, -1.09), "doped Si back gate", fs=6.0)
# thin layers keep short callouts (too thin to hold text)
callout("Au (5 nm)", (XR[1], 2.6, 0.5 * (ZTI + ZAU)), (0.872, 0.600))
callout("Ti (60 nm)", (XR[1], 2.3, 0.5 * (ZTA + ZTI)), (0.872, 0.525))
callout("Ta (10 nm,\nadhesion)", (XR[1], 2.0, 0.5 * (ZG + ZTA)),
        (0.872, 0.430))
callout("monolayer graphene", (1.8, GY[0], 0.035),
        (0.030, 0.068), ha="left")
# supercurrent through the junction ------------------------------------
q0, q1 = p2(4.00, 3.6, 0.34), p2(6.00, 3.6, 0.34)
ov.annotate("", xy=q1, xytext=q0,
            arrowprops=dict(arrowstyle="-|>", lw=1.1, color=C[1],
                            mutation_scale=9))
ov.text(0.5 * (q0[0] + q1[0]) + 0.004, 0.5 * (q0[1] + q1[1]) - 0.055,
        r"$I_s(\varphi,T_e)$", fontsize=6.2, color=C[1], ha="center",
        bbox=dict(fc="white", ec="none", alpha=0.55, pad=0.6))
# quasiparticle exchange: channel <-> contact ---------------------------
e0, e1 = p2(5.55, 5.1, 0.13), p2(6.85, 5.1, 0.60)
ov.annotate("", xy=e1, xytext=e0,
            arrowprops=dict(arrowstyle="<|-|>", lw=0.9, color=C[0],
                            mutation_scale=7))
ov.text(0.5 * (e0[0] + e1[0]) + 0.085, 0.5 * (e0[1] + e1[1]) + 0.080,
        "quasiparticle\nexchange, $\\tau_{\\rm A}$", fontsize=FS,
        color=C[0], ha="center")
# incident microwave photon --------------------------------------------
ph1 = p2(4.65, 4.5, 0.11)
phx, phy = 0.255, 0.870
tt = np.linspace(0, 1, 160)
wob = 0.011 * np.sin(2 * np.pi * 6.5 * tt) * (1 - 0.65 * tt)
dxp, dyp = ph1[0] - phx, ph1[1] - phy
nn = np.hypot(dxp, dyp)
ux, uy = dxp / nn, dyp / nn
px = phx + dxp * tt - uy * wob
py = phy + dyp * tt + ux * wob
ov.plot(px[:-8], py[:-8], color=C[4], lw=0.9)
ov.annotate("", xy=ph1, xytext=(px[-9], py[-9]),
            arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C[4],
                            mutation_scale=7))
ov.text(phx - 0.010, phy + 0.030, "microwave photon, $h\\nu$",
        fontsize=FS, color=C[4], ha="left")
# electron-phonon cooling ----------------------------------------------
c0, c1 = p2(6.55, 0.9, 0.06), p2(6.55, 0.9, -1.30)
ov.annotate("", xy=c1, xytext=c0,
            arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C[2],
                            mutation_scale=7))
ov.text(c0[0] + 0.022, 0.5 * (c0[1] + c1[1]),
        "e-ph cooling\n$\\Sigma A\\,(T_e^3-T_p^3)$", fontsize=FS,
        color=shade(C[2], 0.85), ha="left", va="center")
# dimensions ------------------------------------------------------------
d0, d1 = p2(XL[1], 0.82, 0.02), p2(XR[0], 0.82, 0.02)
ov.annotate("", xy=d1, xytext=d0,
            arrowprops=dict(arrowstyle="<|-|>", lw=0.7, color="0.1",
                            mutation_scale=6))
ov.text(0.5 * (d0[0] + d1[0]) + 0.005, 0.5 * (d0[1] + d1[1]) - 0.052,
        r"$L=0.2\ \mu$m", fontsize=FS, ha="center")
w0, w1 = p2(11.0, YC[0], -1.45), p2(11.0, YC[1], -1.45)
ov.annotate("", xy=w1, xytext=w0,
            arrowprops=dict(arrowstyle="<|-|>", lw=0.7, color="0.1",
                            mutation_scale=6))
ang = np.degrees(np.arctan2(w1[1] - w0[1], w1[0] - w0[0]))
ov.text(0.5 * (w0[0] + w1[0]) + 0.034, 0.5 * (w0[1] + w1[1]) - 0.040,
        r"$W=5.3\ \mu$m", fontsize=FS, ha="center", rotation=ang,
        rotation_mode="anchor")
# readout resonator (top right) ----------------------------------------
r0 = p2(9.15, 5.85, ZAU)
rx, ry = 0.800, 0.915
ov.plot([r0[0], rx - 0.055], [r0[1], ry], color="0.15", lw=0.7)
th = np.linspace(np.pi, 0, 40)
for k in range(4):
    ov.plot(rx - 0.049 + 0.026 * k + 0.013 * np.cos(th),
            ry + 0.016 * np.sin(th), color="0.15", lw=0.7)
ov.plot([rx + 0.055, rx + 0.072], [ry, ry], color="0.15", lw=0.7)
for xc in (rx + 0.072, rx + 0.084):
    ov.plot([xc, xc], [ry - 0.016, ry + 0.016], color="0.15", lw=0.7)
ov.plot([rx + 0.084, rx + 0.100], [ry, ry], color="0.15", lw=0.7)
ov.text(rx + 0.024, ry + 0.048,
        "readout resonator\n$\\nu_r=6$ GHz,  $L_r=2$ nH",
        fontsize=FS, ha="center", va="center")
# back-gate terminal ----------------------------------------------------
g0 = p2(10.0, 1.1, -1.25)
gx_, gy_ = 0.930, 0.075
ov.plot([g0[0], gx_ - 0.014], [g0[1], gy_], color="0.15", lw=0.7)
circ = plt.Circle((gx_, gy_), 0.014, fc="white", ec="0.15", lw=0.7,
                  transform=ov.transAxes)
ov.add_patch(circ)
ov.text(gx_ + 0.022, gy_, r"$V_{\rm BG}$", fontsize=FS, va="center")
ov.text(0.012, 0.012, "not to scale", fontsize=5.0, color="0.45")

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig_device.pdf"))
print("fig_device done")
