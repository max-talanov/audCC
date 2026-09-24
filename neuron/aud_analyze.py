"""
Analysis of auditory-input runs of ctx_thalamus_mpi.py (PLAN-auditory-input.md
design D6). Reads the saved .npz; no NEURON or MPI needed.

    # evoked responses; --sham is the same run without --stim, so each tone
    # window has an exact no-tone counterfactual (the network is deterministic)
    python3 neuron/aud_analyze.py run_tones.npz --sham run_notones.npz --outdir out

    # regression / rank-count check: per-population spike counts of two runs
    python3 neuron/aud_analyze.py --compare legacy.npz nrem.npz

Reports:
  - per-population PSTH around tone onsets, evoked spike counts in 0-50 and
    50-300 ms, paired against the sham run (Wilcoxon signed-rank);
  - TC first-spike latency and the fraction of evoked TC spikes in bursts;
  - TC firing mode for the whole run: burst fraction (Lu et al. 1992: >= 100
    ms silence, then ISI <= 4 ms);
  - cortical state: ISI CV, pairwise spike-count correlation (50 ms bins),
    SO-band (0.5-2 Hz) fraction of the LFP-proxy power, fraction of time in
    cortical population silences >= 100 ms (DOWN states).
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ctx_analyze as A                                  # noqa: E402

POPS = ["tc", "re", "l4e", "l4i", "l23e", "l23i", "l5e", "l5i", "l6e", "l6i"]
CX_E = ["l4e", "l23e", "l5e", "l6e"]
SKIP_MS = 1000.0   # start-up transient excluded from state metrics


def load(path):
    d = np.load(path, allow_pickle=True)
    return {k: d[k] for k in d.files}


def _sel(d, pop):
    lo, hi = d["ranges"].item()[pop]
    m = (d["gids"] >= lo) & (d["gids"] < hi)
    return d["times"][m], d["gids"][m], lo, hi


def burst_mask(t, g, silence=100.0, isi=4.0):
    """Per-spike bool: spike belongs to a Lu et al. (1992) low-threshold burst
    (first spike after >= `silence` ms without spikes, next ISI <= `isi`,
    burst continues while ISI <= `isi`)."""
    order = np.lexsort((t, g))
    ts, gs = t[order], g[order]
    out = np.zeros(len(ts), bool)
    for cell in np.unique(gs):
        idx = np.nonzero(gs == cell)[0]
        x = ts[idx]
        i = 0
        while i < len(x):
            prev_gap = x[i] - x[i - 1] if i > 0 else np.inf
            if prev_gap >= silence and i + 1 < len(x) and x[i + 1] - x[i] <= isi:
                j = i
                while j + 1 < len(x) and x[j + 1] - x[j] <= isi:
                    j += 1
                out[idx[i:j + 1]] = True
                i = j + 1
            else:
                i += 1
    res = np.zeros(len(t), bool)
    res[order] = out
    return res


def psth(t, n_cells, onsets, win=(-300.0, 700.0), bin_ms=5.0):
    edges = np.arange(win[0], win[1] + bin_ms, bin_ms)
    c = np.zeros(len(edges) - 1)
    for t0 in onsets:
        c += np.histogram(t - t0, edges)[0]
    rate = c / (len(onsets) * n_cells * bin_ms / 1000.0)
    return 0.5 * (edges[:-1] + edges[1:]), rate


def trial_counts(t, onsets, a, b):
    return np.array([((t >= t0 + a) & (t < t0 + b)).sum() for t0 in onsets], float)


def evoked_table(d, sham, onsets):
    from scipy.stats import wilcoxon
    rows = []
    for pop in POPS:
        t, _, lo, hi = _sel(d, pop)
        n = hi - lo
        row = {"pop": pop, "n": n}
        for a, b in [(0.0, 50.0), (50.0, 300.0)]:
            c = trial_counts(t, onsets, a, b) / n
            row["stim_%d_%d" % (a, b)] = c.mean()
            if sham is not None:
                ts, _, _, _ = _sel(sham, pop)
                cs = trial_counts(ts, onsets, a, b) / n
                row["sham_%d_%d" % (a, b)] = cs.mean()
                diff = c - cs
                try:
                    p = wilcoxon(c, cs).pvalue if np.any(diff != 0) else 1.0
                except ValueError:
                    p = float("nan")
                row["p_%d_%d" % (a, b)] = p
        rows.append(row)
    return rows


def tc_evoked(d, onsets, win=100.0):
    t, g, _, _ = _sel(d, "tc")
    bm = burst_mask(t, g)
    lat, nb, ne = [], 0, 0
    for t0 in onsets:
        m = (t >= t0) & (t < t0 + win)
        if m.any():
            lat.append(t[m].min() - t0)
        m50 = (t >= t0) & (t < t0 + 50.0)
        ne += m50.sum()
        nb += (m50 & bm).sum()
    return {"lat_median": float(np.median(lat)) if lat else float("nan"),
            "lat_iqr": (float(np.percentile(lat, 25)), float(np.percentile(lat, 75)))
            if lat else (float("nan"),) * 2,
            "trials_with_spike": len(lat), "n_trials": len(onsets),
            "evoked_burst_frac": nb / ne if ne else float("nan")}


def state_metrics(d, rng_seed=0):
    tstop = float(d["tstop"])
    out = {}
    t, g, _, _ = _sel(d, "tc")
    m = t >= SKIP_MS
    bm = burst_mask(t[m], g[m])
    out["tc_rate"] = m.sum() / ((tstop - SKIP_MS) / 1000.0) / (_sel(d, "tc")[3] - _sel(d, "tc")[2])
    out["tc_burst_frac"] = bm.mean() if len(bm) else float("nan")

    # cortical excitatory cells
    ts, gs = [], []
    for pop in CX_E:
        t, g, _, _ = _sel(d, pop)
        m = t >= SKIP_MS
        ts.append(t[m]); gs.append(g[m])
    t, g = np.concatenate(ts), np.concatenate(gs)
    cvs = []
    for cell in np.unique(g):
        x = np.sort(t[g == cell])
        if len(x) >= 5:
            isi = np.diff(x)
            cvs.append(isi.std() / isi.mean())
    out["cx_isi_cv"] = float(np.median(cvs)) if cvs else float("nan")

    edges = np.arange(SKIP_MS, tstop + 50.0, 50.0)
    cells = [c for c in np.unique(g) if (g == c).sum() >= 5]
    rng = np.random.default_rng(rng_seed)
    if len(cells) > 100:
        cells = list(rng.choice(cells, 100, replace=False))
    if len(cells) >= 2:
        X = np.array([np.histogram(t[g == c], edges)[0] for c in cells], float)
        X = X[X.std(axis=1) > 0]
        C = np.corrcoef(X)
        out["cx_pair_corr"] = float(C[np.triu_indices_from(C, 1)].mean())
    else:
        out["cx_pair_corr"] = float("nan")

    counts = np.histogram(t, np.arange(SKIP_MS, tstop + 1.0, 1.0))[0]
    silent, run = 0, 0
    for c in counts:
        if c == 0:
            run += 1
        else:
            if run >= 100:
                silent += run
            run = 0
    if run >= 100:
        silent += run
    out["cx_silence_frac"] = silent / len(counts)

    from scipy.signal import welch
    _, comp, _, bins, _ = A._composite_lfp(d)
    x = comp[bins >= SKIP_MS]
    f, p = welch(x, fs=1000.0, nperseg=min(len(x), 8192))
    band = lambda lo, hi: p[(f >= lo) & (f < hi)].sum()
    out["so_power_frac"] = band(0.5, 2.0) / band(0.5, 40.0)
    return out


def figure(d, sham, onsets, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(POPS) // 2, 2, figsize=(10, 11), sharex=True)
    for ax, pop in zip(axes.T.ravel(), POPS):
        t, _, lo, hi = _sel(d, pop)
        x, r = psth(t, hi - lo, onsets)
        ax.plot(x, r, color=A.RED, lw=1.2, label="tones")
        if sham is not None:
            ts, _, _, _ = _sel(sham, pop)
            x, rs = psth(ts, hi - lo, onsets)
            ax.plot(x, rs, color=A.GREY, lw=1.0, label="sham (no tone)")
        ax.axvspan(0, 50, color=A.BLUE, alpha=0.12)
        ax.set_title(pop.upper(), fontsize=9, loc="left")
        ax.set_ylabel("Hz/cell", fontsize=8)
        ax.tick_params(labelsize=7)
    axes[0, 0].legend(fontsize=7, frameon=False)
    for ax in axes[-1]:
        ax.set_xlabel("time from tone onset (ms)")
    fig.suptitle("Tone-evoked PSTH (%d tones, 50 ms, state=%s); blue = tone"
                 % (len(onsets), str(d.get("state", "?"))), fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


LFP_ROWS = [("l23e", "L2/3"), ("l4e", "L4"), ("l5e", "L5"), ("l6e", "L6"),
            ("tc", "TC"), ("re", "RE")]
LFP_SIGN = {"l5e": -1.0, "l6e": -1.0, "re": -1.0}   # as ctx_analyze._composite_lfp


def _lfp_pair(d, sham, fs=1000.0):
    """Per-population LFP proxy (ctx_analyze._synaptic_lfp, same sign
    convention as _composite_lfp) for the tone run and the sham run, both
    z-scored with the SHAM's mean and SD so their amplitudes are comparable.
    Returns (bins, {pop: tone trace}, {pop: sham trace}); the "ctx" entry is
    the composite cortical LFP (mean of the four layers)."""
    tstop = float(d["tstop"])
    out_t, out_s = {}, {}
    bins = None
    for pop, _ in LFP_ROWS:
        tt, _, _, _ = _sel(d, pop)
        lt, bins = A._synaptic_lfp(tt, tstop, fs=fs)
        if sham is not None:
            ts, _, _, _ = _sel(sham, pop)
            ls, _ = A._synaptic_lfp(ts, tstop, fs=fs)
        else:
            ls = lt
        ref = ls[bins >= SKIP_MS]
        mu, sd = ref.mean(), ref.std() + 1e-9
        sign = LFP_SIGN.get(pop, 1.0)
        out_t[pop] = sign * (lt - mu) / sd
        out_s[pop] = sign * (ls - mu) / sd
    for dct in (out_t, out_s):
        dct["ctx"] = np.mean([dct[p] for p in CX_E], axis=0)
    return bins, out_t, out_s


def lfp_figure(d, sham, onsets, out_png, window=None, pre=300.0, post=1000.0):
    """(a) composite cortical + thalamic LFP proxy over a window of several
    tones, tone run vs sham; (b) tone-triggered average LFP per layer, tone
    vs sham, mean +/- SEM over tones; (c) tone-triggered average of the
    0.5-2 Hz slow-oscillation band of the composite LFP."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fs = 1000.0
    bins, lt, ls = _lfp_pair(d, sham, fs)
    tstop = float(d["tstop"])
    if window is None:
        t0 = onsets[min(3, len(onsets) - 1)] - 1500.0
        window = (max(SKIP_MS, t0), min(tstop, t0 + 12000.0))
    w = (bins >= window[0]) & (bins < window[1])
    tsec = bins[w] / 1000.0

    fig = plt.figure(figsize=(13, 13))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 1.6, 0.8], hspace=0.35,
                          top=0.95)

    ax = fig.add_subplot(gs[0])
    for key, lab, off, col in [("ctx", "cortex (L2/3-L6)", 0.0, "0.2"),
                               ("tc", "TC", -5.0, A.BLUE), ("re", "RE", -10.0, A.RED)]:
        if sham is not None:
            ax.plot(tsec, ls[key][w] + off, color=A.GREY, lw=0.7, alpha=0.8)
        ax.plot(tsec, lt[key][w] + off, color=col, lw=0.9)
        ax.text(tsec[0] - 0.1, off, lab, ha="right", va="center", fontsize=8)
    for t0 in onsets[(onsets >= window[0]) & (onsets < window[1])]:
        ax.axvline(t0 / 1000.0, color=A.BLUE, lw=1.0, alpha=0.35)
    ax.set_yticks([])
    ax.set_xlabel("time (s)")
    ax.set_title("(a) LFP proxy, %.1f-%.1f s: tone run (colour) vs sham (grey); "
                 "vertical lines = tone onsets" % (window[0] / 1000, window[1] / 1000),
                 fontsize=10, loc="left")

    lags = np.arange(-pre, post, 1000.0 / fs)
    idx0 = [int(round(t0 * fs / 1000.0)) for t0 in onsets]
    ok = [i for i in idx0 if i - pre >= 0 and i + post < len(bins)]

    def epochs(x):
        return np.array([x[i - int(pre):i + int(post)] for i in ok])

    sub = gs[1].subgridspec(2, 3, hspace=0.45, wspace=0.25)
    for k, (pop, lab) in enumerate(LFP_ROWS):
        ax = fig.add_subplot(sub[k // 3, k % 3])
        for tr, col, name in [(ls, A.GREY, "sham"), (lt, A.RED, "tones")]:
            if tr is ls and sham is None:
                continue
            E = epochs(tr[pop])
            m, se = E.mean(0), E.std(0) / np.sqrt(len(E))
            ax.fill_between(lags, m - se, m + se, color=col, alpha=0.25, lw=0)
            ax.plot(lags, m, color=col, lw=1.1, label=name)
        ax.axvspan(0, 50, color=A.BLUE, alpha=0.12)
        ax.axhline(0, color="k", lw=0.4)
        ax.set_title(lab, fontsize=9, loc="left")
        ax.tick_params(labelsize=7)
        if k % 3 == 0:
            ax.set_ylabel("z (sham SD)", fontsize=8)
        if k >= 3:
            ax.set_xlabel("time from tone onset (ms)", fontsize=8)
        if k == 0:
            ax.legend(fontsize=7, frameon=False)
    fig.axes[1].text(0.0, 1.28, "(b) Tone-triggered average LFP proxy per population, "
                     "mean ± SEM over %d tones (blue = tone)" % len(ok),
                     transform=fig.axes[1].transAxes, fontsize=10)

    ax = fig.add_subplot(gs[2])
    for tr, col, name in [(ls, A.GREY, "sham"), (lt, "#8e44ad", "tones")]:
        if tr is ls and sham is None:
            continue
        so = A._bandpass(tr["ctx"], fs, 0.5, 2.0)
        E = epochs(so)
        m, se = E.mean(0), E.std(0) / np.sqrt(len(E))
        ax.fill_between(lags, m - se, m + se, color=col, alpha=0.25, lw=0)
        ax.plot(lags, m, color=col, lw=1.2, label=name)
    ax.axvspan(0, 50, color=A.BLUE, alpha=0.12)
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlabel("time from tone onset (ms)")
    ax.set_ylabel("z")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("(c) Tone-triggered 0.5-2 Hz (slow-oscillation) band of the cortical "
                 "LFP proxy -- an evoked slow wave would show here", fontsize=10, loc="left")

    fig.suptitle("Tone-evoked LFP proxy (state=%s, %d tones of %s dB, 50 ms). LFP proxy "
                 "= spikes x PSP kernel (rise 2 ms, decay 100 ms), not a volume-conductor LFP"
                 % (str(d.get("state", "?")), len(onsets),
                    "/".join("%g" % x for x in np.unique(d["stim_level"]))),
                 fontsize=10, y=0.985)
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)


def reconstruction_figure(d, sham, onsets, out_png, window=(20000.0, 30000.0)):
    """ctx_analyze's literature-style format (Fernandez & Luthi 2020 Fig. 1A
    style): SO-band + 10-15 Hz band of the composite cortical LFP proxy, the
    10-15 Hz band below it, detected 10-15 Hz events shaded grey
    (ctx_analyze._detect_spindles). One pair of panels for the tone run and
    one for the sham, on identical y-scales, with tone onsets in blue."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fs = 1000.0
    cases = [("NREM + tones (60 dB, 50 ms)", d, "#b03a2e")]
    if sham is not None:
        cases.append(("NREM, no tones (sham)", sham, "0.35"))
    comp = []
    for label, npz, color in cases:
        _, composite, _, bins, _ = A._composite_lfp(npz, fs=fs)
        so = A._bandpass(composite, fs, 0.5, 2.0)
        spin = A._bandpass(composite, fs, 10.0, 15.0)
        starts, ends = A._detect_spindles(spin, bins)
        comp.append((label, color, bins, so + spin, spin, starts, ends))
    w0 = lambda b: (b >= window[0]) & (b < window[1])
    raw_max = max(np.abs(c[3][w0(c[2])]).max() for c in comp)
    spin_max = max(np.abs(c[4][w0(c[2])]).max() for c in comp)
    fig, axes = plt.subplots(2 * len(comp), 1, figsize=(14, 2.3 * 2 * len(comp) + 0.6),
                             sharex=True)
    tone_s = onsets[(onsets >= window[0]) & (onsets < window[1])] / 1000.0
    for i, (label, color, bins, rec, spin, starts, ends) in enumerate(comp):
        w = w0(bins)
        tsec = bins[w] / 1000.0
        ax_raw, ax_spin = axes[2 * i], axes[2 * i + 1]
        ax_raw.plot(tsec, rec[w], color=color, lw=0.8)
        ax_raw.set_ylim(-raw_max * 1.1, raw_max * 1.1)
        ax_raw.set_ylabel("raw\n(SO+spindle)", fontsize=8)
        ax_raw.set_title(label, fontsize=10, loc="left", fontweight="bold")
        ax_spin.plot(tsec, spin[w], color=color, lw=0.8)
        ax_spin.set_ylim(-spin_max * 1.1, spin_max * 1.1)
        ax_spin.set_ylabel("10-15 Hz", fontsize=8)
        for s_, e_ in zip(starts, ends):
            a, b = bins[s_] / 1000.0, bins[e_] / 1000.0
            if b < window[0] / 1000.0 or a > window[1] / 1000.0:
                continue
            for ax in (ax_raw, ax_spin):
                ax.axvspan(a, b, color="0.5", alpha=0.18, lw=0)
        for ax in (ax_raw, ax_spin):
            for t in tone_s:
                ax.axvline(t, color=A.BLUE, lw=1.2, alpha=0.8, ls="--" if i else "-")
    axes[-1].set_xlabel("time (s)")
    axes[-1].set_xlim(window[0] / 1000.0, window[1] / 1000.0)
    fig.suptitle("Cortical LFP proxy, SO band (0.5-2 Hz) + spindle band (10-15 Hz); "
                 "grey = detected 10-15 Hz events, blue = tone onsets "
                 "(dashed in the sham: the same times, no tone)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def compare(a_path, b_path):
    a, b = load(a_path), load(b_path)
    print("%-6s %10s %10s %8s" % ("pop", Path(a_path).stem, Path(b_path).stem, "diff %"))
    worst = 0.0
    for pop in POPS:
        na, nb = len(_sel(a, pop)[0]), len(_sel(b, pop)[0])
        d = 100.0 * (nb - na) / na if na else (0.0 if nb == 0 else float("inf"))
        worst = max(worst, abs(d))
        print("%-6s %10d %10d %+8.2f" % (pop, na, nb, d))
    # spike rasters as (time, gid) rows sorted by time then gid: independent
    # of the order ranks were gathered in
    def raster(d):
        order = np.lexsort((d["gids"], np.round(d["times"], 6)))
        return np.column_stack([np.round(d["times"], 6)[order], d["gids"][order]])
    ka, kb = raster(a), raster(b)
    same = ka.shape == kb.shape and np.array_equal(ka, kb)
    n = min(len(ka), len(kb))
    neq = np.nonzero(np.any(ka[:n] != kb[:n], axis=1))[0]
    first = ka[neq[0], 0] if len(neq) else (None if len(ka) == len(kb) else ka[n - 1, 0])
    print("identical spike rasters: %s; first divergence at t = %s ms; "
          "largest per-population count difference %.2f%%"
          % (same, "none" if first is None else "%.1f" % first, worst))
    return worst, same


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz", nargs="?")
    ap.add_argument("--sham", default="")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--compare", nargs=2, metavar=("A", "B"))
    ap.add_argument("--recon-windows", type=float, nargs="*", default=[20000.0],
                    help="start times (ms) of 10 s windows for the "
                         "literature-style SO+spindle reconstruction figure")
    a = ap.parse_args()
    if a.compare:
        compare(*a.compare)
        return
    d = load(a.npz)
    sham = load(a.sham) if a.sham else None
    tag = Path(a.npz).stem
    if "stim_t" in d:
        tstop = float(d["tstop"])
        onsets = d["stim_t"][d["stim_t"] + 300.0 < tstop]
        print("== %s: %d tones, state=%s ==" % (tag, len(onsets), d.get("state")))
        rows = evoked_table(d, sham, onsets)
        hdr = "%-5s %6s | %-24s | %-24s" % ("pop", "cells", "spikes/cell 0-50 ms",
                                             "spikes/cell 50-300 ms")
        print(hdr)
        print("%-5s %6s | %-24s | %-24s" % ("", "", "tone  sham   p", "tone  sham   p"))
        for r in rows:
            cols = []
            for w in ("0_50", "50_300"):
                if sham is not None:
                    cols.append("%5.2f %5.2f %7.1e" % (r["stim_" + w], r["sham_" + w], r["p_" + w]))
                else:
                    cols.append("%5.2f" % r["stim_" + w])
            print("%-5s %6d | %-24s | %-24s" % (r["pop"], r["n"], cols[0], cols[1]))
        te = tc_evoked(d, onsets)
        print("TC first-spike latency: median %.1f ms (IQR %.1f-%.1f), %d/%d tones "
              "with a TC spike in 0-100 ms; evoked TC spikes in bursts: %.0f%%"
              % (te["lat_median"], te["lat_iqr"][0], te["lat_iqr"][1],
                 te["trials_with_spike"], te["n_trials"], 100 * te["evoked_burst_frac"]))
        os.makedirs(a.outdir, exist_ok=True)
        png = os.path.join(a.outdir, "%s_psth.png" % tag)
        figure(d, sham, onsets, png)
        print("figure:", png)
        png = os.path.join(a.outdir, "%s_lfp.png" % tag)
        lfp_figure(d, sham, onsets, png)
        print("figure:", png)
        for w0 in a.recon_windows:
            png = os.path.join(a.outdir, "%s_reconstructed_%d-%ds.png"
                               % (tag, w0 / 1000, (w0 + 10000) / 1000))
            reconstruction_figure(d, sham, onsets, png, window=(w0, w0 + 10000.0))
            print("figure:", png)
    sm = state_metrics(d)
    print("state: TC %.2f Hz/cell, TC burst fraction %.2f | cortex ISI CV %.2f, "
          "pair corr %.3f, silence>=100ms %.2f, SO-band power frac %.2f"
          % (sm["tc_rate"], sm["tc_burst_frac"], sm["cx_isi_cv"],
             sm["cx_pair_corr"], sm["cx_silence_frac"], sm["so_power_frac"]))


if __name__ == "__main__":
    main()
