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
    sm = state_metrics(d)
    print("state: TC %.2f Hz/cell, TC burst fraction %.2f | cortex ISI CV %.2f, "
          "pair corr %.3f, silence>=100ms %.2f, SO-band power frac %.2f"
          % (sm["tc_rate"], sm["tc_burst_frac"], sm["cx_isi_cv"],
             sm["cx_pair_corr"], sm["cx_silence_frac"], sm["so_power_frac"]))


if __name__ == "__main__":
    main()
