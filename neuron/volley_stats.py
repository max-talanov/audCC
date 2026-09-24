"""
Are the "spindles" in a ctx_thalamus_mpi.py run real multi-cycle spindles, or
single thalamic volleys that the 10-15 Hz band-pass turns into a short ripple?

For each run (.npz spike dump) this reports:
  - thalamic (TC+RE) volleys: count, rate, peak, inter-volley interval + CV
  - spindle-like trains: fraction of big-peak intervals < 200 ms
  - follow-up activity 50-250 ms after a volley, relative to the volley peak
  - cortical activity between volleys, and which layers carry it
  - filter-ringing control: spindle-band envelope width around volleys in the
    composite LFP proxy vs. a pure impulse train at the same times through
    the same filter (equal widths => the "spindle" is the filter's response)
and optionally plots thalamus/cortex population rates for a short window.

Usage:
    python3 neuron/volley_stats.py \
        "before=res/2026-08-31/ctx_nrn_45171023.npz" \
        "Option 2=res/2026-09-07/ctx_nrn_45453403.npz" \
        --plot out/mn5_population_rates_before_vs_option2.png
"""

import argparse
import os
import sys

import numpy as np
from scipy.signal import find_peaks, hilbert

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ctx_analyze as A                                    # noqa: E402

THAL = ("tc", "re")
CTX_E = ("l23e", "l4e", "l5e", "l6e")


def _pop_rate(t, g, ranges, pops, tstop, smooth_ms=3.0, skip_ms=1000.0):
    """Population rate in Hz/cell, 1 ms bins; first skip_ms zeroed."""
    m = np.zeros_like(t, bool)
    for k in pops:
        m |= (g >= ranges[k][0]) & (g < ranges[k][1])
    n = sum(ranges[k][1] - ranges[k][0] for k in pops)
    r, b = A._rate(t[m], tstop, 1.0, smooth_ms)
    r = r / n * 1000.0
    r[b < skip_ms] = 0.0
    return r, b


def volley_stats(npz, volley_hz=100.0, min_gap_ms=200.0, skip_ms=1000.0):
    t, g, R, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    r, b = _pop_rate(t, g, R, THAL, tstop, skip_ms=skip_ms)
    pk, pr = find_peaks(r, height=volley_hz, distance=int(min_gap_ms))
    heights = pr["peak_heights"]
    ivi = np.diff(b[pk])
    follow = np.array([r[i + 50:i + 250].max() / r[i] for i in pk[:-1]])
    pk2, _ = find_peaks(r, height=0.3 * np.median(heights), distance=40)
    trains = float(np.mean(np.diff(b[pk2]) < 200.0))

    # which layers fire far (>150 ms) from any thalamic volley
    vt = b[pk]
    between = {}
    for k in CTX_E + ("l5i",):
        tk = t[(g >= R[k][0]) & (g < R[k][1]) & (t > skip_ms)]
        idx = np.searchsorted(vt, tk)
        prev = vt[np.clip(idx - 1, 0, len(vt) - 1)]
        nxt = vt[np.clip(idx, 0, len(vt) - 1)]
        far = (np.abs(tk - prev) > 150) & (np.abs(nxt - tk) > 150)
        between[k] = far.sum() / (R[k][1] - R[k][0]) / ((tstop - skip_ms) / 1000.0)

    # filter-ringing control on the composite LFP proxy
    _, comp, _, bins, _ = A._composite_lfp(npz)
    spin = A._bandpass(comp, 1000.0, 10.0, 15.0)
    imp = np.zeros_like(comp)
    imp[pk] = 1.0
    spin_imp = A._bandpass(imp, 1000.0, 10.0, 15.0)
    lags = np.arange(-300, 600)
    inner = pk[(pk + lags[0] >= 0) & (pk + lags[-1] < len(comp))]

    def halfwidth(x):
        env = np.abs(hilbert(x))
        prof = np.mean([env[i + lags] for i in inner], axis=0)
        a = lags[prof > prof.max() / 2]
        return float(a.max() - a.min())

    return dict(n=len(pk), rate=len(pk) / ((tstop - skip_ms) / 1000.0),
                peak=float(np.median(heights)), ivi_median=float(np.median(ivi)),
                ivi_p5=float(np.percentile(ivi, 5)), ivi_p95=float(np.percentile(ivi, 95)),
                ivi_cv=float(ivi.std() / ivi.mean()), trains=trains,
                follow=float(np.median(follow)),
                active=float(np.mean(r[b >= skip_ms] > 20.0)),
                between_hz=between,
                halfwidth_model=halfwidth(spin), halfwidth_impulse=halfwidth(spin_imp))


def window_table(npz, windows_ms):
    """Per-window stats for one run, to check the dynamics are stationary:
    thalamic volleys, their mean interval + CV, detected 'spindle' events
    (ctx_analyze._detect_spindles, threshold from the whole run), SO-band
    RMS of the composite LFP proxy, RE mean burst size, and L5 RS / IB
    firing rates (IB = first half of the L5E gid range)."""
    t, g, R, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    _, comp, _, bins, _ = A._composite_lfp(npz)
    so = A._bandpass(comp, 1000.0, 0.5, 2.0)
    starts, _ = A._detect_spindles(A._bandpass(comp, 1000.0, 10.0, 15.0), bins)
    r, b = _pop_rate(t, g, R, THAL, tstop)
    vt = b[find_peaks(r, height=100.0, distance=200)[0]]
    lo, hi = R["l5e"]
    nib = int(round((hi - lo) * 0.5))
    rows = []
    for w0, w1 in windows_ms:
        sec = (w1 - w0) / 1000.0
        mv = (vt >= w0) & (vt < w1)
        iv = np.diff(vt[mv])
        tw = (t >= w0) & (t < w1)
        re = A._re_burst_stats(t[tw], g[tw], *R["re"])
        rows.append(dict(
            window=(w0, w1), volleys=int(mv.sum()), volley_rate=mv.sum() / sec,
            ivi_mean=float(iv.mean()) if len(iv) else float("nan"),
            ivi_cv=float(iv.std() / iv.mean()) if len(iv) > 1 else float("nan"),
            spindles=int(((bins[starts] >= w0) & (bins[starts] < w1)).sum()),
            so_rms=float(np.sqrt(np.mean(so[(bins >= w0) & (bins < w1)] ** 2))),
            re_burst=re["mean_burst_size"],
            l5_rs=((g >= lo + nib) & (g < hi) & tw).sum() / (hi - lo - nib) / sec,
            l5_ib=((g >= lo) & (g < lo + nib) & tw).sum() / nib / sec))
    return rows


def plot_rates(cases, out_png, window=(20000, 23000)):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2 * len(cases), 1, figsize=(14, 2.5 * 2 * len(cases)), sharex=True)
    for i, (label, npz) in enumerate(cases):
        t, g, R, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
        for j, (pops, name) in enumerate([(THAL, "thalamus TC+RE"), (CTX_E, "cortex E (L2/3-L6)")]):
            r, b = _pop_rate(t, g, R, pops, tstop)
            w = (b >= window[0]) & (b < window[1])
            ax = axes[2 * i + j]
            ax.plot(b[w] / 1000.0, r[w], lw=0.6, color=A.BLUE if j == 0 else A.RED)
            ax.set_ylabel("Hz/cell", fontsize=8)
            ax.set_title(f"{label}: {name}", fontsize=9, loc="left")
    axes[-1].set_xlabel("time (s)")
    fig.suptitle("Population firing rate (1 ms bins, 3 ms smoothing) -- single thalamic volleys, "
                 "no multi-cycle spindle trains", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    print(f"Saved population-rate figure to {out_png}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help='"label=path.npz" pairs')
    ap.add_argument("--plot", default=None, help="optional output PNG for the rate plot")
    ap.add_argument("--window-start", type=float, default=20000.0)
    ap.add_argument("--window-len", type=float, default=3000.0)
    ap.add_argument("--windows", default="",
                    help='per-window table instead of the full-run report, e.g. '
                         '"20-40,60-80,100-120,160-180" (seconds)')
    a = ap.parse_args(argv)

    cases = []
    for pair in a.runs:
        label, path = pair.split("=", 1)
        cases.append((label, np.load(path, allow_pickle=True)))
    if a.windows:
        wins = [tuple(float(x) * 1000.0 for x in w.split("-")) for w in a.windows.split(",")]
        for label, npz in cases:
            print(f"\n== {label}")
            print("  window      | volleys (rate)  | mean interval (CV) | 'spindles' | SO RMS "
                  "| RE burst | L5 RS / IB Hz/cell")
            for row in window_table(npz, wins):
                w0, w1 = row["window"]
                print(f"  {w0/1000:4.0f}-{w1/1000:4.0f} s | {row['volleys']:3d} ({row['volley_rate']:.2f}/s) "
                      f"| {row['ivi_mean']:4.0f} ms ({row['ivi_cv']:.2f})     | {row['spindles']:3d}        "
                      f"| {row['so_rms']:.3f}  | {row['re_burst']:.2f}     | {row['l5_rs']:.1f} / {row['l5_ib']:.1f}")
        return 0
    for label, npz in cases:
        s = volley_stats(npz)
        print(f"\n== {label}")
        print(f"  thalamic volleys: n={s['n']} ({s['rate']:.2f}/s), peak {s['peak']:.0f} Hz/cell")
        print(f"  inter-volley interval: median {s['ivi_median']:.0f} ms "
              f"(5-95% {s['ivi_p5']:.0f}-{s['ivi_p95']:.0f}), CV {s['ivi_cv']:.2f}")
        print(f"  spindle-like trains (big-peak intervals < 200 ms): {100 * s['trains']:.1f}%")
        print(f"  strongest activity 50-250 ms after a volley: {100 * s['follow']:.0f}% of the volley peak")
        print(f"  time thalamus > 20 Hz/cell: {100 * s['active']:.0f}%")
        print("  firing >150 ms from any volley (Hz/cell): "
              + ", ".join(f"{k} {v:.2f}" for k, v in s["between_hz"].items()))
        print(f"  spindle-band envelope half-width around volleys: model {s['halfwidth_model']:.0f} ms "
              f"vs pure impulse {s['halfwidth_impulse']:.0f} ms")
    if a.plot:
        plot_rates(cases, a.plot, window=(a.window_start, a.window_start + a.window_len))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
