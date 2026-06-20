"""
Draw the auditory thalamo-cortical loop architecture as a schematic figure.

This renders the wiring actually built in ``tc_network.connect_network`` /
``attach_sleep_drives``: the feed-forward thalamocortical input, the
intracortical flow, and -- the additions that turn cc's feed-forward column
into a closed sleeping loop -- the corticothalamic feedback (L6 -> MGB/nRT) and
the reciprocal reticulo-thalamic inhibition (nRT <-> TCR) that forms the ~13 Hz
spindle resonator. The 1 Hz / 13 Hz sleep drives and the auditory background are
shown as inputs.

    python3 tc_architecture.py --outdir out

Saves ``out/tc_architecture.png``. Pure matplotlib; no NEST required.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


# ---- palette ---------------------------------------------------------------
EXC = "#c0392b"      # excitatory projections
INH = "#2471a3"      # inhibitory projections
LOOP = "#2f8f6f"     # corticothalamic feedback (loop-closing)
DRIVE = "#7d3c98"    # injected sleep drives / background
CORTEX_FC = "#fdf2e9"
CORTEX_EC = "#e59866"
THAL_FC = "#eaf2f8"
THAL_EC = "#5499c7"
INHBOX_FC = "#eaf2f8"


def _box(ax, xy, w, h, label, fc, ec, fontsize=10, sub=None, lw=1.6):
    x, y = xy
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2))
    ax.text(x, y + (0.10 if sub else 0), label, ha="center", va="center",
            fontsize=fontsize, fontweight="bold", zorder=3)
    if sub:
        ax.text(x, y - 0.22, sub, ha="center", va="center",
                fontsize=fontsize - 2.5, color="0.30", zorder=3)
    return {"c": (x, y), "w": w, "h": h}


def _edge(box, side, frac=0.0):
    x, y = box["c"]
    w, h = box["w"], box["h"]
    if side == "top":
        return (x + frac * w / 2, y + h / 2)
    if side == "bottom":
        return (x + frac * w / 2, y - h / 2)
    if side == "left":
        return (x - w / 2, y + frac * h / 2)
    if side == "right":
        return (x + w / 2, y + frac * h / 2)
    return (x, y)


def _arrow(ax, p0, p1, color, kind="exc", lw=2.0, rad=0.0, ls="-", z=4):
    # excitatory -> filled arrowhead; inhibitory -> bracket (flat-bar) head
    style = "-|>" if kind == "exc" else "-["
    mut = 14 if kind == "exc" else 7
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle=style, mutation_scale=mut,
        connectionstyle=f"arc3,rad={rad}", linewidth=lw, linestyle=ls,
        color=color, zorder=z, shrinkA=2, shrinkB=2))


def build_figure():
    fig, ax = plt.subplots(figsize=(10.5, 12))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 13)
    ax.axis("off")

    cx = 4.4          # cortical column centre
    iw, ih = 2.9, 0.95  # E-box size
    ix = cx + 2.45    # inhibitory interneuron column x

    # ---- cortical excitatory layers (top -> down) ----
    L23 = _box(ax, (cx, 11.4), iw, ih, "L2/3", CORTEX_FC, CORTEX_EC, sub="RS + FRB")
    L4 = _box(ax, (cx, 9.8), iw, ih, "L4", CORTEX_FC, CORTEX_EC, sub="spiny stellate")
    L5 = _box(ax, (cx, 8.2), iw, ih, "L5", CORTEX_FC, CORTEX_EC, sub="TuftRS + TuftIB")
    L6 = _box(ax, (cx, 6.6), iw, ih, "L6", CORTEX_FC, CORTEX_EC, sub="cortico-thalamic")

    # ---- cortical inhibitory interneuron pools ----
    I23 = _box(ax, (ix, 11.4), 1.5, 0.8, "I", INHBOX_FC, INH, fontsize=9,
               sub="Bask/LTS/Axax")
    I4 = _box(ax, (ix, 9.8), 1.5, 0.8, "I", INHBOX_FC, INH, fontsize=9, sub="L4 I")
    I56 = _box(ax, (ix, 7.4), 1.5, 0.8, "I", INHBOX_FC, INH, fontsize=9,
               sub="L5/6 shared")

    # ---- thalamus (MGB) ----
    TCR = _box(ax, (cx - 0.9, 3.4), 2.6, 1.05, "TCR", THAL_FC, THAL_EC,
               sub="MGB relay (E)")
    nRT = _box(ax, (cx + 2.4, 3.4), 2.4, 1.05, "nRT", INHBOX_FC, INH,
               sub="reticular (I)")
    ax.text(cx + 0.6, 4.35, "Thalamus (MGB)", ha="center", fontsize=11,
            fontstyle="italic", color="0.35")
    ax.add_patch(FancyBboxPatch((cx - 2.45, 2.75), 6.1, 1.55,
                 boxstyle="round,pad=0.02,rounding_size=0.1",
                 linewidth=1.2, edgecolor="0.7", facecolor="none",
                 linestyle=(0, (4, 3)), zorder=1))

    # ---- intracortical excitatory flow (L4 -> L2/3 -> L5 -> L6) ----
    _arrow(ax, _edge(L4, "top", -0.3), _edge(L23, "bottom", -0.3), EXC)
    _arrow(ax, _edge(L23, "bottom", 0.3), _edge(L5, "top", -0.15), EXC, rad=-0.32)
    _arrow(ax, _edge(L4, "bottom", 0.25), _edge(L5, "top", 0.25), EXC, rad=-0.15)
    _arrow(ax, _edge(L5, "bottom", -0.25), _edge(L6, "top", -0.25), EXC)
    # L6 <-> L5 recurrent
    _arrow(ax, _edge(L6, "top", 0.3), _edge(L5, "bottom", 0.3), EXC, rad=-0.3, lw=1.6)

    # ---- intracortical inhibition (I -| E, per layer) ----
    _arrow(ax, _edge(I23, "left"), _edge(L23, "right"), INH, kind="inh", lw=1.6)
    _arrow(ax, _edge(I4, "left"), _edge(L4, "right"), INH, kind="inh", lw=1.6)
    _arrow(ax, _edge(I56, "left", 0.4), _edge(L5, "right", -0.3), INH, kind="inh",
           lw=1.6, rad=0.15)
    _arrow(ax, _edge(I56, "left", -0.4), _edge(L6, "right", 0.3), INH, kind="inh",
           lw=1.6, rad=-0.15)
    # E -> I (drive onto interneurons), light
    _arrow(ax, _edge(L23, "right", 0.5), _edge(I23, "left", 0.6), EXC, lw=1.0,
           rad=0.2)
    _arrow(ax, _edge(L4, "right", 0.5), _edge(I4, "left", 0.6), EXC, lw=1.0, rad=0.2)

    # ---- thalamocortical input: TCR -> L4 (dense, thick) ----
    _arrow(ax, _edge(TCR, "top", -0.2), _edge(L4, "bottom", -0.6), EXC, lw=3.0)
    ax.text(cx - 1.95, 6.1, "thalamo-\ncortical", fontsize=8.5, color=EXC,
            ha="center", va="center", fontweight="bold")

    # ---- corticothalamic feedback: L6 -> TCR and L6 -> nRT (closes the loop) ----
    _arrow(ax, _edge(L6, "bottom", -0.55), _edge(TCR, "top", 0.3), LOOP, lw=2.6,
           rad=0.12)
    _arrow(ax, _edge(L6, "bottom", 0.7), _edge(nRT, "top", -0.1), LOOP, lw=2.6,
           rad=-0.18)
    ax.text(cx + 2.15, 5.55, "cortico-thalamic\nfeedback", fontsize=8.5,
            color=LOOP, ha="center", va="center", fontweight="bold")

    # ---- spindle resonator: nRT -| TCR, TCR -> nRT, nRT self-inhibition ----
    _arrow(ax, _edge(nRT, "left", 0.4), _edge(TCR, "right", 0.4), INH, kind="inh",
           lw=2.4, rad=-0.25)
    _arrow(ax, _edge(TCR, "right", -0.4), _edge(nRT, "left", -0.4), EXC,
           lw=2.0, rad=-0.25)
    # nRT self-inhibition loop
    _arrow(ax, _edge(nRT, "bottom", 0.3), _edge(nRT, "bottom", -0.3), INH,
           kind="inh", lw=1.6, rad=-1.6)
    ax.text(cx + 0.75, 2.35, "nRT <-> TCR\nspindle resonator (~13 Hz)", fontsize=8,
            color="0.25", ha="center", va="center", fontweight="bold")

    # ---- drives / background inputs ----
    slow = _box(ax, (8.7, 10.6), 2.0, 0.95, "1 Hz", "#f5eef8", DRIVE,
                fontsize=10, sub="slow osc (AC)")
    spin = _box(ax, (8.7, 3.4), 2.0, 0.95, "13 Hz", "#f5eef8", DRIVE,
                fontsize=10, sub="spindle (AC)")
    aud = _box(ax, (0.95, 3.4), 1.7, 0.95, "MGB in", "#f5eef8", DRIVE,
               fontsize=9, sub="auditory\nPoisson")

    # slow 1 Hz onto cortical pyramidal pools
    _arrow(ax, _edge(slow, "left"), (_edge(L23, "right", 0.9)[0] + 0.1,
           _edge(L23, "right", 0.9)[1]), DRIVE, lw=1.8, ls=(0, (5, 2)), rad=0.0)
    _arrow(ax, _edge(slow, "bottom"), _edge(L5, "top", 0.9), DRIVE, lw=1.5,
           ls=(0, (5, 2)), rad=-0.3)
    # 13 Hz + 1 Hz gating onto thalamus
    _arrow(ax, _edge(spin, "left"), _edge(nRT, "right"), DRIVE, lw=1.8,
           ls=(0, (5, 2)))
    # auditory background onto TCR
    _arrow(ax, _edge(aud, "right"), _edge(TCR, "left"), DRIVE, lw=1.8,
           ls=(0, (5, 2)))

    # ---- legend ----
    lx, ly = 0.55, 12.4
    leg = [
        (EXC, "-|>", "-", "excitatory"),
        (INH, "-[", "-", "inhibitory"),
        (LOOP, "-|>", "-", "cortico-thalamic loop"),
        (DRIVE, "-|>", (0, (5, 2)), "sleep drive / input"),
    ]
    for i, (c, st, ls, lab) in enumerate(leg):
        yy = ly - i * 0.42
        ax.add_patch(FancyArrowPatch((lx, yy), (lx + 0.7, yy), arrowstyle=st,
                     mutation_scale=12, linewidth=2.0, linestyle=ls, color=c,
                     zorder=5))
        ax.text(lx + 0.85, yy, lab, fontsize=9, va="center")

    ax.set_title("Auditory thalamo-cortical loop architecture (NEST sleep model)",
                 fontsize=13, fontweight="bold", pad=12)
    fig.tight_layout()
    return fig


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=str, default="out")
    ap.add_argument("--name", type=str, default="tc_architecture.png")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig = build_figure()
    out = outdir / args.name
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved architecture figure to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
