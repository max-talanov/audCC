"""
Stage 2 report (PLAN-spindels.md): slow oscillation, spike-based spindles
and their nesting, for ctx_thalamus_mpi.py runs (.npz). No NEURON needed.

  - SO: UP-state onsets = peaks of the L5E + L6E population rate (20 ms
    smoothing, >= 300 ms apart, above mean + 2 SD); rate and interval CV
  - spindles: Stage 0 RE volley trains (volley_stats.py); counted as
    spindle-like with >= 6 cycles
  - nesting: a train "starts in an UP state" if its first volley is within
    -50 .. +300 ms of an UP-state onset; fraction of UP states that carry a
    train with >= 2 / >= 6 cycles

    python3 neuron/stage2_report.py "label=run.npz" ... [--skip 2000]
"""

import argparse
import os
import sys

import numpy as np
from scipy.signal import find_peaks

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ctx_analyze as A                                  # noqa: E402
import volley_stats as V                                 # noqa: E402


def up_onsets(t, g, R, tstop, skip):
    m = np.zeros_like(t, bool)
    n = 0
    for k in ("l5e", "l6e"):
        lo, hi = R[k]
        m |= (g >= lo) & (g < hi)
        n += hi - lo
    r, b = A._rate(t[m], tstop, 1.0, 20.0)
    r = r / max(n, 1) * 1000.0
    r[b < skip] = 0.0
    thr = r[b >= skip].mean() + 2 * r[b >= skip].std()
    pk, _ = find_peaks(r, height=thr, distance=300)
    return b[pk]


def report(label, npz, skip=2000.0):
    t, g, R, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    span = (tstop - skip) / 1000.0
    ups = up_onsets(t, g, R, tstop, skip)
    iui = np.diff(ups)
    trains = [tr for tr in V.volley_trains(V.re_volleys(t, g, R, tstop, skip_ms=skip))]
    cyc = np.array([tr["cycles"] for tr in trains]) if trains else np.zeros(0, int)
    sp = [tr for tr in trains if tr["cycles"] >= 6]
    multi = [tr for tr in trains if tr["cycles"] >= 2]

    def in_up(tr):
        return bool(len(ups)) and np.any((tr["t0"] - ups >= -50) & (tr["t0"] - ups <= 300))

    def carried(k):
        return np.mean([any(u - 50 <= tr["t0"] <= u + 300 and tr["cycles"] >= k for tr in trains)
                        for u in ups]) if len(ups) else float("nan")

    print("\n== %s" % label)
    print("  SO: %d UP states (%.2f Hz), interval median %.0f ms, CV %.2f"
          % (len(ups), len(ups) / span, np.median(iui) if len(iui) else np.nan,
             iui.std() / iui.mean() if len(iui) > 1 else np.nan))
    if not trains:
        print("  no RE volleys")
        return
    hist = " ".join("%d:%d" % (k, (cyc == k).sum()) for k in range(1, min(cyc.max(), 15) + 1)
                    if (cyc == k).sum())
    print("  RE trains: %d (%.1f/min); cycles median %.0f, max %d; histogram %s"
          % (len(trains), len(trains) / span * 60, np.median(cyc), cyc.max(), hist))
    if multi:
        print("  multi-cycle trains: freq median %.1f Hz, duration median %.0f ms"
              % (np.median([tr["freq"] for tr in multi]), np.median([tr["duration"] for tr in multi])))
    print("  spindle-like (>= 6 cycles): %d (%.1f/min)%s"
          % (len(sp), len(sp) / span * 60,
             ", %.1f-%.1f Hz, %.0f-%.0f ms" % (min(x["freq"] for x in sp), max(x["freq"] for x in sp),
                                              min(x["duration"] for x in sp), max(x["duration"] for x in sp))
             if sp else ""))
    print("  nesting: trains starting in an UP state %.0f%% (>= 6 cycles: %s); UP states carrying a "
          ">= 2-cycle train %.0f%%, a >= 6-cycle train %.0f%%"
          % (100 * np.mean([in_up(tr) for tr in trains]),
             "%.0f%%" % (100 * np.mean([in_up(tr) for tr in sp])) if sp else "-",
             100 * carried(2), 100 * carried(6)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help='"label=path.npz" pairs')
    ap.add_argument("--skip", type=float, default=2000.0)
    a = ap.parse_args()
    for pair in a.runs:
        label, path = pair.split("=", 1)
        report(label, np.load(path, allow_pickle=True), a.skip)


if __name__ == "__main__":
    main()
