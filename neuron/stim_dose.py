"""
E0 dose-response (PLAN-stimulation-experiments.md WP1): response to a single
stimulus as a function of level, stim vs sham, from --open-loop runs.

Per level: TC and L4E spikes/cell in the first 50 ms, P(UP onset within 300 ms),
P(spindle within 300 ms). Sham fake tones at the same scheduled times and
levels are the control. Wilson 95% CIs; Fisher / Mann-Whitney p.

    python3 neuron/stim_dose.py --pair stim.npz=sham.npz [--pair ...] \\
        --cache rows.pkl --out-png figures/e0_dose_response.png --out-txt table.txt
"""

import argparse
import os
import pickle

import numpy as np

import aud_phase as AP
import cl_analyze as CLA
import stim_phase as SPH


def collect(pairs, skip=2000.0, up_thr=2.0):
    rows = {"stim": [], "sham": []}
    for k, pair in enumerate(pairs):
        for arm, path in zip(("stim", "sham"), pair):
            z = np.load(path, allow_pickle=True)
            t, info = SPH.tones_of_run(z)
            for r, i in zip(AP.measure(z, t, skip, up_thr), info):
                r.update(i)
                r["seed"] = k
                rows[arm].append(r)
    return rows


def table(rows, levels):
    f = lambda rs, k: np.array([r[k] for r in rs], float)
    out = [f"{'level dB':>8} | {'n stim/sham':>11} | {'TC sp/cell 50ms':>20} | {'L4E sp/cell 50ms':>20} | "
           f"{'P(UP<300ms)':>24} | {'P(spindle<300ms)':>24}", "-" * 120]
    for lv in list(levels) + ["all"]:
        a = [r for r in rows["stim"] if lv == "all" or r["level"] == lv]
        b = [r for r in rows["sham"] if lv == "all" or r["level"] == lv]
        if not a or not b:
            continue
        cells = []
        for k in ("tc", "l4e"):
            cells.append(f"{f(a,k).mean():5.2f} vs {f(b,k).mean():4.2f} (p={CLA._p_cont(f(a,k), f(b,k)):.3f})")
        for k in ("evoked_up", "evoked_sp"):
            cells.append(f"{f(a,k).mean():4.0%} vs {f(b,k).mean():4.0%} (p={CLA._p_bool(f(a,k), f(b,k)):.3f})")
        out.append(f"{str(lv):>8} | {len(a):>4}/{len(b):<6} | " + " | ".join(f"{c:>20}" for c in cells))
    return out


def figure(rows, levels, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [("tc", "TC spikes/cell, 0-50 ms", False), ("l4e", "L4 E spikes/cell, 0-50 ms", False),
              ("evoked_up", "P(UP onset within 300 ms)", True),
              ("evoked_sp", "P(spindle within 300 ms)", True)]
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6))
    x = np.arange(len(levels))
    for ax, (key, title, binary) in zip(axs, panels):
        for arm, col, off in (("sham", "tab:gray", -0.1), ("stim", "tab:red", 0.1)):
            m, lo, hi = [], [], []
            for lv in levels:
                v = np.array([r[key] for r in rows[arm] if r["level"] == lv], float)
                m.append(v.mean())
                if binary:
                    a, b = CLA._wilson(v.sum(), len(v))
                else:
                    se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
                    a, b = v.mean() - 1.96 * se, v.mean() + 1.96 * se
                lo.append(a), hi.append(b)
            m, lo, hi = map(np.array, (m, lo, hi))
            ax.errorbar(x + off, m, yerr=[np.clip(m - lo, 0, None), np.clip(hi - m, 0, None)],
                        fmt="o-", color=col, capsize=3, lw=1.4, label=arm)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{int(l)}" for l in levels])
        ax.set_xlabel("stimulus level (dB)")
        ax.set_title(title, fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
        if binary:
            ax.set_ylim(-0.02, 1.05)
    axs[0].legend(frameon=False, fontsize=8)
    n = len(rows["stim"])
    fig.suptitle(f"E0 dose-response, one 50 ms stimulus, {n} stim / {len(rows['sham'])} sham tones; "
                 "bars = 95% CI", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", action="append", required=True)
    ap.add_argument("--cache", default="")
    ap.add_argument("--out-png", default="")
    ap.add_argument("--out-txt", default="")
    a = ap.parse_args(argv)
    if a.cache and os.path.exists(a.cache):
        rows = pickle.load(open(a.cache, "rb"))
    else:
        rows = collect([p.split("=", 1) for p in a.pair])
        if a.cache:
            pickle.dump(rows, open(a.cache, "wb"))
    levels = sorted({r["level"] for r in rows["stim"] + rows["sham"]})
    txt = "\n".join(table(rows, levels))
    print(txt)
    if a.out_txt:
        open(a.out_txt, "w").write(txt + "\n")
    if a.out_png:
        figure(rows, levels, a.out_png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
