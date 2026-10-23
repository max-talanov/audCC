"""
Per-cell-type PSTHs locked to the tone (M2 in PLAN-stimulation-experiments.md):
population firing rate (Hz per cell, 20 ms bins) of TC, RE and the E/I
populations of L4, L2/3, L5, L6, stim vs sham, split by causal SO phase at the
tone (stim_phase.py) and by spindle presence. Mean +/- SEM over tones.

    python3 neuron/stim_psth.py --pair stim.npz=sham.npz [--pair ...] --out-prefix figures/stim_psth
"""

import argparse

import numpy as np

import cl_analyze as CLA
import stim_phase as SPH

POPS = [("tc", "TC"), ("re", "RE"), ("l4e", "L4 E"), ("l4i", "L4 I"), ("l23e", "L2/3 E"),
        ("l23i", "L2/3 I"), ("l5e", "L5 E"), ("l5i", "L5 I"), ("l6e", "L6 E"), ("l6i", "L6 I")]
PRE, POST, BIN = 500.0, 1500.0, 20.0
EDGES = np.arange(-PRE, POST + BIN, BIN)


def arm_psth(run, rows):
    """{pop: tones x bins rate in Hz/cell} for the tones in rows."""
    z = run.z
    t, g, R = z["times"], z["gids"], z["ranges"].item()
    out = {}
    for name, _ in POPS:
        lo, hi = R[name]
        ts = np.sort(t[(g >= lo) & (g < hi)])
        rate = np.empty((len(rows), len(EDGES) - 1))
        for k, r in enumerate(rows):
            c = np.searchsorted(ts, r["t"] + EDGES)
            rate[k] = np.diff(c) / (hi - lo) / (BIN / 1000.0)
        out[name] = rate
    return out


def collect(pairs, skip=2000.0):
    S = {"stim": [], "sham": []}
    for sp, sh in pairs:
        for arm, path in (("stim", sp), ("sham", sh)):
            run = SPH.Run(path)
            tones, info = SPH.tones_of_run(run.z)
            rows = [r for r in SPH.tone_rows(run, tones, skip, info=info)
                    if r["t"] - PRE >= 0 and r["t"] + POST <= run.T]
            S[arm].append((rows, arm_psth(run, rows)))
    return S


def pool(S, sel):
    acc = {n: [] for n, _ in POPS}
    for rows, P in S:
        idx = [j for j, r in enumerate(rows) if sel(r)]
        for n, _ in POPS:
            if idx:
                acc[n].append(P[n][idx])
    return {n: (np.concatenate(v) if v else np.empty((0, len(EDGES) - 1))) for n, v in acc.items()}


def figure(S, cols, out_png, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tt = EDGES[:-1] + BIN / 2
    fig, axs = plt.subplots(len(POPS), len(cols), figsize=(3.0 * len(cols), 1.7 * len(POPS)),
                            sharex=True, squeeze=False)
    for j, (lab, sel) in enumerate(cols):
        P = {arm: pool(S[arm], sel) for arm in ("stim", "sham")}
        for i, (name, nm) in enumerate(POPS):
            ax = axs[i, j]
            for arm, col in (("sham", "tab:gray"), ("stim", "tab:red")):
                v = P[arm][name]
                if len(v) < 3:
                    continue
                m, se = v.mean(0), v.std(0, ddof=1) / np.sqrt(len(v))
                ax.plot(tt, m, color=col, lw=1.1, label=f"{arm} (n={len(v)})" if i == 0 else None)
                ax.fill_between(tt, m - se, m + se, color=col, alpha=0.25, lw=0)
            ax.axvline(0, color="k", lw=0.6, ls=":")
            ax.set_ylim(bottom=0)
            ax.spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax.set_ylabel(f"{nm}\nHz/cell", fontsize=7)
            if i == 0:
                ax.set_title(lab, fontsize=9)
                ax.legend(frameon=False, fontsize=6, loc="upper right")
    for ax in axs[-1]:
        ax.set_xlabel("time from tone (ms)", fontsize=8)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", action="append", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--split", default="phase,spindle",
                    help="comma list of column splits: phase, spindle, level (open-loop runs)")
    a = ap.parse_args(argv)
    S = collect([p.split("=", 1) for p in a.pair])
    ttl = ("Tone-locked population rate per cell type (20 ms bins), mean +/- SEM over tones\n"
           "stim = real tone, sham = fake tone at the same point of the closed loop")
    sp = a.split.split(",")
    if "phase" in sp:
        cols = [("all tones", lambda r: True)] + [
            (f"SO {n}", (lambda b: lambda r: r["so_bin"] == b)(b)) for b, n in enumerate(SPH.SO_NAMES)]
        figure(S, cols, a.out_prefix + "_by_so_phase.png", ttl + "\ncolumns: SO phase at the tone (causal estimate)")
    if "spindle" in sp:
        cols = [("no spindle at tone", lambda r: not r["sp_present"]),
                ("spindle present at tone", lambda r: r["sp_present"])]
        figure(S, cols, a.out_prefix + "_by_spindle.png", ttl + "\ncolumns: spindle at the tone")
    if "level" in sp:
        levels = sorted({r["level"] for arm in S.values() for rows, _ in arm for r in rows})
        cols = [("all tones", lambda r: True)] + [
            (f"{lv:g} dB", (lambda v: lambda r: r["level"] == v)(lv)) for lv in levels]
        figure(S, cols, a.out_prefix + "_by_level.png", ttl + "\ncolumns: stimulus level")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
