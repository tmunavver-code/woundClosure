import sys, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection
sys.path.insert(0, "/home/claude/phase")
from vxport import C_HEX, unwrap

c = C_HEX
df = pd.read_csv("/home/claude/phase/out/sweep_results.csv")
SOFT_BG, HEX_BG = "#e3f0f9", "#cfcfcf"
INK, MUTED = "#1f1f1f", "#6b6b6b"
HEXPT, SOFTPT = "#3a3a3a", "#1f6fb5"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK})

fig = plt.figure(figsize=(17.5, 7.4))
ax = fig.add_axes([0.045, 0.20, 0.47, 0.72])
ax2 = fig.add_axes([0.745, 0.20, 0.245, 0.72])
Gm, Lmin, Lmax = 0.18, -1.4, 0.6

# --- theory regions (Farhadifar et al. 2007) --------------------------------
Gc = np.linspace(0, 2 / c**2, 300)
LR = 8 / c * ((1 - 0.5 * Gc * c**2) / 3) ** 1.5          # pressure-free collapse boundary
ax.add_patch(Polygon([[Lmin, 0], [0, 0], [-2 * c * Gm, Gm], [Lmin, Gm]], color=SOFT_BG, lw=0, zorder=0))
hexpoly = [[0, 0], [-2 * c * Gm, Gm], [0, Gm], [0, 2 / c**2]] + list(zip(LR[::-1], Gc[::-1]))
ax.add_patch(Polygon(hexpoly, color=HEX_BG, lw=0, zorder=0))
ax.add_patch(Polygon([[0, Gm], [0, 2 / c**2]] + list(zip(LR[::-1], Gc[::-1])) + [[Lmax, 0], [Lmax, Gm]],
                     facecolor="white", edgecolor="#9a9a9a", hatch="//", lw=0, zorder=0))
ax.plot([0, -2 * c * Gm], [0, Gm], color=INK, lw=1.8, zorder=2)
ax.plot(list(LR) + [0, 0], list(Gc) + [2 / c**2, Gm], color=INK, lw=1.3, ls="--", zorder=2)

# --- simulated points -------------------------------------------------------
on = df[(df.Lbar <= Lmax + 1e-9) & ~df.set.str.startswith("point_") & (df.set != "inset_soft")]
h = on[on.phase == "hexagonal"]; s = on[on.phase == "soft"]
ax.scatter(h.Lbar, h.Gbar, marker="s", s=22, color=HEXPT, edgecolor="white", lw=0.5, zorder=3,
           label=f"simulated: hexagonal ({len(h)})")
ax.scatter(s.Lbar, s.Gbar, marker="o", s=24, color=SOFTPT, edgecolor="white", lw=0.5, zorder=3,
           label=f"simulated: soft ({len(s)})")
pts = df[df.set.str.startswith("point_")]
for _, r in pts.iterrows():
    name = r.set.split("_")[1]
    ax.scatter(r.Lbar, r.Gbar, s=150, color="#8cc63f", edgecolor=INK, lw=1.0, zorder=4)
    ax.scatter(r.Lbar, r.Gbar, s=12, color=INK, zorder=5)
    right = name == "I"
    ax.annotate(f"{name}  ({r.phase})", (r.Lbar, r.Gbar), xytext=(9 if right else -9, 7), textcoords="offset points",
                ha="left" if right else "right", fontsize=10.5, color=INK, zorder=6,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))

ax.text(-1.32, 0.052, "Soft Network", fontsize=13, weight="bold", color=INK, bbox=dict(fc=SOFT_BG, ec="none", pad=2))
ax.text(-0.52, 0.163, "Hexagonal Network", fontsize=13, weight="bold", color=INK, bbox=dict(fc=HEX_BG, ec="none", pad=2))
ax.text(0.27, 0.118, "collapse\n(pressure-free box;\nnot tested: model\nbox is fixed)", fontsize=8.5,
        color=MUTED, bbox=dict(fc="white", ec="none", pad=1.5))
ax.set_xlim(Lmin, Lmax); ax.set_ylim(0, Gm)
ax.set_xlabel(r"normalized line tension  $\bar\Lambda=\Lambda/K A_0^{3/2} = -2p_0/(k_a r_{stiff})$")
ax.set_ylabel(r"normalized contractility  $\bar\Gamma=\Gamma/K A_0 = 1/(k_a r_{stiff})$")
ax.set_title("Ground-state phase diagram reproduced with the 2D_VX_model energy and forces",
             fontsize=12, color=INK, loc="left")
ax.plot([], [], color=INK, lw=1.8, label=r"theory: $p_0=\sqrt{8\sqrt{3}}=3.7224$")
ax.plot([], [], color=INK, lw=1.3, ls="--", label="theory: collapse (pressure-free)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2, fontsize=9.5, frameon=False)


# --- insets: relaxed networks from the simulation ----------------------------
def draw_net(axi, f):
    d = np.load(f)
    X = unwrap(d["r"], d["conn"], float(d["Lx"]), float(d["Ly"]))
    polys = [x for x in X]
    axi.add_collection(PolyCollection(polys, facecolor="white", edgecolor=INK, lw=0.9))
    axi.set_xlim(0, float(d["Lx"])); axi.set_ylim(0, float(d["Ly"]))
    axi.set_aspect("equal"); axi.set_xticks([]); axi.set_yticks([])
    for sp in axi.spines.values(): sp.set_color(INK)
    return d

ins1 = fig.add_axes([0.545, 0.575, 0.135, 0.31])
d1 = draw_net(ins1, "/home/claude/phase/out/cfg_II.npz")
qII = df.loc[df.set == "point_II", "q_mean"].iloc[0]
qIII = df.loc[df.set == "point_III", "q_mean"].iloc[0]
ins1.set_title(f"relaxed network, point II\nhexagonal, q = {qII:.4f}", fontsize=9, color=INK, pad=3)
ins2 = fig.add_axes([0.545, 0.20, 0.135, 0.31])
d2 = draw_net(ins2, "/home/claude/phase/out/cfg_III.npz")
ins2.set_title(f"relaxed network, point III\nsoft, q = {qIII:.4f}  (p$_0$ = 4.25)", fontsize=9, color=INK, pad=3)
ax2.text(3.775, 3.728, "all six $\\bar\\Gamma$ rows lie on\none curve (overlapping)", fontsize=9.5, color=MUTED)

# --- right panel: going along / across the line -------------------------------
line = df[df.set == "line"].copy()
Gs = sorted(line.Gbar.unique())
cmap = plt.get_cmap("Blues")
for i, G in enumerate(Gs):
    sub = line[np.isclose(line.Gbar, G)].sort_values("p0")
    col = cmap(0.35 + 0.65 * i / max(1, len(Gs) - 1))
    ax2.plot(sub.p0, sub.q_mean, "-o", color=col, ms=5, lw=1.6, label=rf"$\bar\Gamma$ = {G:.2f}")
pp = np.linspace(3.53, 3.9, 200)
ax2.plot(pp, np.maximum(pp, c), color=INK, lw=1, ls=":", label=r"theory: $q=\max(p_0,\,3.7224)$")
ax2.axvline(c, color=INK, lw=1.0, ls="--")
ax2.text(c + 0.004, 3.86, "p$_0^*$ = 3.7224", fontsize=9.5, color=INK)
ax2.set_xlabel(r"preferred shape index $p_0$  (crossing the line at fixed $\bar\Gamma$)")
ax2.set_ylabel(r"relaxed mean shape index  $\langle P/\sqrt{A}\rangle$")
ax2.set_title("Along the line: all six rows switch at p$_0$ = 3.7224", fontsize=12, color=INK, loc="left")
ax2.grid(color="#e6e6e6", lw=0.8); ax2.set_axisbelow(True)
ax2.legend(fontsize=9, loc="upper left", frameon=True, facecolor="white", edgecolor="#bbbbbb")
for sp in ["top", "right"]: ax2.spines[sp].set_visible(False)

fig.savefig("/mnt/user-data/outputs/phase_diagram_validation.png", dpi=170, facecolor="white")
print("saved")
