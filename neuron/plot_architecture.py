"""
Architecture diagrams of the thalamocortical model as PNG (figures/):

  figures/architecture_cortex.png          cortical column: 4 layers, E / I
  figures/architecture_thalamocortical.png TC <-> RE loop, L4 / L6, inputs
  figures/architecture_traub.png           Traub et al. 2005 cell types vs this model

Neuron and synapse counts are those of the full-scale configuration
(SCALE=1.65, SPINDLE=stage2 + --cx-noise-rate 200 + K_Na in RE), counted by
building ParallelCorticoThalamicNet and counting connections per projection
(2026-10-12). Orange = excitatory, blue = inhibitory, grey = external input.

    python3 neuron/plot_architecture.py [--outdir figures]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch     # noqa: E402
from matplotlib.path import Path                                     # noqa: E402

EXC = ("#FAC775", "#854F0B", "#633806")    # fill, edge, text
INH = ("#B5D4F4", "#185FA5", "#0C447C")
EXT = ("#D3D1C7", "#5F5E5A", "#2C2C2A")
A_EXC, A_INH, A_EXT = "#BA7517", "#185FA5", "#888780"
TXT = "#3d3d3a"


def setup(h, w=680):
    fig = plt.figure(figsize=(w / 100 * 1.4, h / 100 * 1.4), dpi=150)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, col, lines):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=8",
                                fc=col[0], ec=col[1], lw=0.8))
    n = len(lines)
    step = h / (n + 1)
    for i, s in enumerate(lines):
        ax.text(x + w / 2, y + step * (i + 1), s, ha="center", va="center",
                fontsize=10 if i == 0 else 8.5, color=col[2],
                fontweight="bold" if i == 0 else "normal")


def arrow(ax, pts, color, dashed=False):
    codes = [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1)
    ax.add_patch(FancyArrowPatch(path=Path(pts, codes), arrowstyle="-|>", mutation_scale=10,
                                 color=color, lw=1.1, ls="--" if dashed else "-",
                                 shrinkA=0, shrinkB=0))


def curve(ax, p0, c1, c2, p1, color):
    path = Path([p0, c1, c2, p1], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=10, color=color,
                                 lw=1.1, shrinkA=0, shrinkB=0))


def label(ax, x, y, s, ha="left"):
    ax.text(x, y, s, ha=ha, va="baseline", fontsize=8.5, color=TXT)


def legend(ax, y, items):
    x = items[0][0]
    for x, col, s in items:
        ax.add_patch(FancyBboxPatch((x, y - 10), 12, 12, boxstyle="round,pad=0,rounding_size=2",
                                    fc=col[0], ec=col[1], lw=0.6))
        label(ax, x + 18, y, s)


def cortex(path):
    fig, ax = setup(470)
    ax.text(20, 16, "Cortical column (scale 1.65): 4,594 neurons", fontsize=11,
            color=TXT, fontweight="bold", va="center")
    rows = [("L2/3", 1056, 264, ["RS pyramidal"]),
            ("L4", 1056, 264, ["RS pyramidal"]),
            ("L5", 874, 214, ["437 IB + 437 RS", "IB gap junctions 3,481"]),
            ("L6", 693, 173, ["RS, corticothalamic"])]
    ei = [("26.4k", "105.6k"), ("26.4k", "105.6k"), ("21.4k", "87.4k"), ("17.3k", "69.3k")]
    for i, (name, ne, ni, sub) in enumerate(rows):
        y = 40 + 100 * i
        box(ax, 190, y, 180, 60, EXC, ["%s E · %d" % (name, ne)] + sub)
        box(ax, 480, y, 140, 60, INH, ["%s I · %d" % (name, ni), "FS interneurons"])
        arrow(ax, [(370, y + 22), (480, y + 22)], A_EXC)
        arrow(ax, [(480, y + 42), (370, y + 42)], A_INH)
        label(ax, 425, y + 15, ei[i][0], "center")
        label(ax, 425, y + 56, ei[i][1], "center")
    arrow(ax, [(240, 140), (240, 100)], A_EXC)
    label(ax, 248, 125, "105.6k")
    arrow(ax, [(320, 200), (320, 240)], A_EXC)
    label(ax, 328, 224, "87.4k")
    arrow(ax, [(240, 300), (240, 340)], A_EXC)
    label(ax, 248, 324, "69.3k")
    curve(ax, (330, 300), (330, 330), (360, 330), (360, 300), A_EXC)
    label(ax, 368, 330, "NMDA recurrent 87.4k")
    arrow(ax, [(190, 55), (110, 55), (110, 285), (190, 285)], A_EXC)
    label(ax, 104, 175, "87.4k", "right")
    arrow(ax, [(190, 260), (150, 260), (150, 80), (190, 80)], A_EXC)
    label(ax, 156, 125, "105.6k")
    legend(ax, 446, [(190, EXC, "Excitatory"), (290, INH, "Inhibitory")])
    label(ax, 390, 446, "Arrow numbers = synapses")
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def loop(path):
    fig, ax = setup(500)
    ax.text(20, 16, "Thalamocortical loop (scale 1.65): 5,031 neurons, 1.16 M synapses",
            fontsize=11, color=TXT, fontweight="bold", va="center")
    for y, h, name in ((30, 130, "Cortex"), (230, 140, "Thalamus")):
        ax.add_patch(FancyBboxPatch((40, y), 600, h, boxstyle="round,pad=0,rounding_size=16",
                                    fc="none", ec="#b4b2a9", lw=0.8, ls="--"))
        ax.text(56, y + 20, name, fontsize=10, color=TXT, fontweight="bold")
    box(ax, 80, 70, 180, 60, EXC, ["L4 E · 1056", "thalamic input layer"])
    box(ax, 420, 70, 180, 60, EXC, ["L6 E · 693", "corticothalamic output"])
    arrow(ax, [(260, 100), (420, 100)], A_EXC, dashed=True)
    label(ax, 340, 90, "via L2/3 and L5", "center")
    box(ax, 80, 270, 180, 70, EXC, ["TC · 346", "relay, I_T + Ca-dep. I_h", "ring, local wiring"])
    box(ax, 420, 270, 180, 70, INH, ["RE (TRN) · 91", "I_T, SK2, K_Na", "lateral 358, gaps 715"])
    arrow(ax, [(260, 290), (420, 290)], A_EXC)
    label(ax, 340, 282, "3,458 (38 per RE)", "center")
    arrow(ax, [(420, 320), (260, 320)], A_INH)
    label(ax, 340, 338, "3,460 (10 per TC)", "center")
    arrow(ax, [(140, 270), (140, 130)], A_EXC)
    label(ax, 148, 206, "105.6k")
    arrow(ax, [(560, 130), (560, 270)], A_EXC)
    label(ax, 568, 206, "9.1k")
    arrow(ax, [(460, 130), (224, 270)], A_EXC)
    label(ax, 330, 182, "34.6k", "right")
    box(ax, 80, 400, 180, 44, EXT, ["IC fibres · 346"])
    arrow(ax, [(170, 400), (170, 340)], A_EXT)
    label(ax, 178, 378, "2,768 (8 per TC)")
    box(ax, 420, 400, 180, 44, EXT, ["Tone stimulus"])
    arrow(ax, [(420, 422), (260, 422)], A_EXT)
    legend(ax, 476, [(80, EXC, "Excitatory"), (180, INH, "Inhibitory"),
                     (280, EXT, "External input (only with --stim)")])
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def traub(path):
    """Traub et al. (2005) cell inventory (J Neurophysiol 93: 2194, Fig. 1
    legend: 14 types, 3,560 neurons) against this model's populations.
    Stage numbers are those of PLAN-architecture.md."""
    # (group, Traub type, n in Traub, E/I, this model, status, note)
    rows = [
        ("L2/3", "SupRS  ·  RS pyramid", 1000, "E", "L2/3 E · 1056", "present", ""),
        ("", "SupFRB  ·  fast rhythmic bursting pyramid", 50, "E", "not modelled", "missing", "Stage 3"),
        ("", "SupBask  ·  basket", 90, "I", "L2/3 I · 264", "present", "one FS type"),
        ("", "SupAxAx  ·  axo-axonic", 90, "I", "not modelled", "missing", "out of scope"),
        ("", "SupLTS  ·  low-threshold spiking", 90, "I", "not modelled", "missing", "Stage 2"),
        ("L4", "SpinStell  ·  spiny stellate", 240, "E", "L4 E · 1056", "partial", "pyramidal stand-in"),
        ("L5", "TuftIB  ·  intrinsically bursting pyramid", 800, "E", "L5 IB · 437", "present", ""),
        ("", "TuftRS  ·  regular spiking pyramid", 200, "E", "L5 RS · 437", "present", ""),
        ("L6", "NontuftRS  ·  pyramid, projects to thalamus", 500, "E", "L6 E · 693", "present", ""),
        ("Deep", "DeepBask  ·  basket", 100, "I", "L5 I · 214 + L6 I · 173", "present", ""),
        ("", "DeepAxAx  ·  axo-axonic", 100, "I", "not modelled", "missing", "out of scope"),
        ("", "DeepLTS  ·  low-threshold spiking", 100, "I", "not modelled", "missing", "Stage 2"),
        ("Thal.", "TCR  ·  thalamocortical relay", 100, "E", "TC · 346", "present", ""),
        ("", "nRT  ·  reticular", 100, "I", "RE · 91", "present", ""),
    ]
    W, H, top, step, bh = 700, 650, 92, 33, 27
    fig, ax = setup(H, W)
    ax.text(20, 16, "Traub et al. 2005 cell types against this model: 9 of 14 present, 5 missing",
            fontsize=11, color=TXT, fontweight="bold", va="center")
    ax.text(20, 34, "Traub: 3,560 neurons (paper Fig. 1).  This model: 5,031 at scale 1.65.",
            fontsize=8.5, color=TXT, va="center")
    xl, wl, xr, wr = 95, 300, 455, 225
    ax.text(xl, top - 14, "Traub et al. 2005  (cells)", fontsize=9, color=TXT, fontweight="bold")
    ax.text(xr, top - 14, "This model", fontsize=9, color=TXT, fontweight="bold")
    for i, (grp, name, n, ei, ours, status, note) in enumerate(rows):
        y = top + i * step
        col = EXC if ei == "E" else INH
        if grp:
            if i:
                ax.plot([20, W - 20], [y - 3, y - 3], color="#d3d1c7", lw=0.6)
            ax.text(22, y + bh / 2, grp, fontsize=9, color=TXT, fontweight="bold", va="center")
        ax.add_patch(FancyBboxPatch((xl, y), wl, bh, boxstyle="round,pad=0,rounding_size=6",
                                    fc=col[0], ec=col[1], lw=0.8))
        ax.text(xl + 8, y + bh / 2, name, fontsize=8.3, color=col[2], va="center")
        ax.text(xl + wl - 8, y + bh / 2, str(n), fontsize=8.5, color=col[2], va="center",
                ha="right", fontweight="bold")
        if status == "missing":
            ax.add_patch(FancyBboxPatch((xr, y), wr, bh, boxstyle="round,pad=0,rounding_size=6",
                                        fc="white", ec="#888780", lw=0.9, ls="--"))
            ax.text(xr + 8, y + bh / 2, ours, fontsize=8.3, color="#888780", va="center")
            ax.text(xr + wr - 8, y + bh / 2, note, fontsize=8, color="#5f5e5a", va="center",
                    ha="right", style="italic")
        else:
            ec = "#c77d00" if status == "partial" else col[1]
            ax.add_patch(FancyBboxPatch((xr, y), wr, bh, boxstyle="round,pad=0,rounding_size=6",
                                        fc=col[0], ec=ec, lw=1.6 if status == "partial" else 0.8,
                                        ls="--" if status == "partial" else "-"))
            ax.text(xr + 8, y + bh / 2, ours, fontsize=8.3, color=col[2], va="center",
                    fontweight="bold")
            if note:
                ax.text(xr + wr - 8, y + bh / 2, note, fontsize=8, color=col[2], va="center",
                        ha="right", style="italic")
            ax.plot([xl + wl, xr], [y + bh / 2] * 2, color=col[1], lw=0.8)
    y0 = top + len(rows) * step + 6
    notes = ["Also in this model but not in Traub: L4 I · 264 (Traub has no layer-4 interneurons;",
             "layer-4 inhibition there comes from the deep interneurons).",
             "Connection classes still missing (PLAN-architecture.md): E to E recurrence outside L5,",
             "any I to I synapse, TC to interneurons and to L5/L6, GABA_B."]
    for k, t in enumerate(notes):
        label(ax, 22, y0 + 14 * k + 8, t)
    legend(ax, H - 16, [(22, EXC, "Excitatory"), (112, INH, "Inhibitory")])
    ax.add_patch(FancyBboxPatch((212, H - 26), 12, 12, boxstyle="round,pad=0,rounding_size=2",
                                fc="#fff", ec="#c77d00", lw=1.4, ls="--"))
    label(ax, 230, H - 16, "stand-in")
    ax.add_patch(FancyBboxPatch((292, H - 26), 12, 12, boxstyle="round,pad=0,rounding_size=2",
                                fc="#fff", ec="#888780", lw=0.9, ls="--"))
    label(ax, 310, H - 16, "missing (planned stage)")
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     "..", "figures"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    for name, fn in (("architecture_cortex.png", cortex),
                     ("architecture_thalamocortical.png", loop),
                     ("architecture_traub.png", traub)):
        p = os.path.join(a.outdir, name)
        fn(p)
        print("Saved", os.path.normpath(p))


if __name__ == "__main__":
    main()
