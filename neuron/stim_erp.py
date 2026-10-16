"""
Event-related field potentials for the stimulation experiments (M1 in
PLAN-stimulation-experiments.md): the LFP proxy of each cortical layer and of
the thalamus, locked to the tone, stim vs sham, split by SO phase at the tone
and by spindle presence.

LFP proxy per population = spikes convolved with a PSP kernel
(ctx_analyze._synaptic_lfp), sign as in ctx_analyze._composite_lfp (deep layers
and RE flipped). Both arms of a seed are normalised with the SHAM run's mean
and SD, so amplitudes are comparable between arms. Each ERP is baseline-
corrected on [-200, 0) ms. The thalamic traces are spike-based proxies, not
measured fields. Tone phase/spindle classes come from stim_phase.py (causal).

    python3 neuron/stim_erp.py --pair stim.npz=sham.npz [--pair ...] --out-prefix figures/stim_erp
"""

import argparse

import numpy as np

import ctx_analyze as A
import stim_phase as SPH
import cl_analyze as CLA

CH = [("l23e", "L2/3"), ("l4e", "L4"), ("l5e", "L5"), ("l6e", "L6"), ("tc", "TC (thalamus)"),
      ("re", "RE (thalamus)")]
PRE, POST, BASE = 500, 1500, 200


def channels(npz, ref=None):
    """name -> LFP proxy at 1 kHz; normalised with ref's (mean, sd) if given."""
    t, g, R, T = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    out, stats = {}, {}
    for name, _ in CH:
        lo, hi = R[name]
        m = (g >= lo) & (g < hi)
        lfp, _ = A._synaptic_lfp(t[m], T)
        sign = -1.0 if name in ("l5e", "l6e", "re") else 1.0
        mu, sd = ref[name] if ref else (lfp.mean(), lfp.std() + 1e-9)
        stats[name] = (mu, sd)
        out[name] = sign * (lfp - mu) / sd
    return out, stats


def erps(lfp, rows):
    """(kept rows, {channel: tones x time array}), baseline-corrected."""
    n = len(lfp[CH[0][0]])
    kept = [r for r in rows if int(r["t"]) - PRE >= 0 and int(r["t"]) + POST <= n]
    out = {}
    for name, _ in CH:
        x = lfp[name]
        out[name] = np.array([x[int(r["t"]) - PRE:int(r["t"]) + POST]
                              - x[int(r["t"]) - BASE:int(r["t"])].mean() for r in kept])
    return kept, out


def collect(pairs, skip=2000.0):
    S = {"stim": [], "sham": []}
    for sp, sh in pairs:
        runs = {"sham": SPH.Run(sh), "stim": SPH.Run(sp)}
        lf_sham, st = channels(runs["sham"].z)
        lf = {"sham": lf_sham, "stim": channels(runs["stim"].z, st)[0]}
        for arm in ("stim", "sham"):
            tones, _ = SPH.tones_of_run(runs[arm].z)
            rows = SPH.tone_rows(runs[arm], tones, skip)
            S[arm].append(erps(lf[arm], rows))
    return S


def pool(S, arm, sel):
    """concatenate tones over seeds for rows satisfying sel(row)."""
    acc = {n: [] for n, _ in CH}
    for rows, E in S:
        for n, _ in CH:
            acc[n].extend(E[n][j] for j, r in enumerate(rows) if sel(r))
    return {n: np.array(v) for n, v in acc.items()}


def figure(S, cols, out_png, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tt = np.arange(-PRE, POST)
    fig, axs = plt.subplots(len(CH), len(cols), figsize=(3.0 * len(cols), 1.9 * len(CH)),
                            sharex=True, squeeze=False)
    for j, (lab, sel) in enumerate(cols):
        P = {arm: pool(S[arm], arm, sel) for arm in ("stim", "sham")}
        for i, (name, nm) in enumerate(CH):
            ax = axs[i, j]
            for arm, col in (("sham", "tab:gray"), ("stim", "tab:red")):
                v = P[arm][name]
                if len(v) < 3:
                    continue
                m, se = v.mean(0), v.std(0, ddof=1) / np.sqrt(len(v))
                ax.plot(tt, m, color=col, lw=1.2, label=f"{arm} (n={len(v)})" if i == 0 else None)
                ax.fill_between(tt, m - se, m + se, color=col, alpha=0.25, lw=0)
            ax.axvline(0, color="k", lw=0.6, ls=":")
            ax.spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax.set_ylabel(nm, fontsize=8)
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
    a = ap.parse_args(argv)
    S = collect([p.split("=", 1) for p in a.pair])
    ttl = ("ERP of the LFP proxy per layer and thalamus, tone-locked\n"
           "mean +/- SEM over tones; z-scored to the sham run, baseline [-200,0) ms")
    cols = [("all tones", lambda r: True)] + [
        (f"SO {n}", (lambda b: lambda r: r["so_bin"] == b)(b)) for b, n in enumerate(SPH.SO_NAMES)]
    figure(S, cols, a.out_prefix + "_by_so_phase.png", ttl + "\ncolumns: SO phase at the tone (causal estimate)")
    cols = [("no spindle at tone", lambda r: not r["sp_present"]),
            ("spindle present at tone", lambda r: r["sp_present"])]
    figure(S, cols, a.out_prefix + "_by_spindle.png", ttl + "\ncolumns: spindle at the tone")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
