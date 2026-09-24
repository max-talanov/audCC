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

Spike-based spindle measure (PLAN-spindels.md Stage 0, --trains): RE
population volleys are detected from spikes, consecutive volleys 60-150 ms
apart are chained into one train, and each train is reported with its cycle
count, duration, within-train frequency and TC / RE participation per cycle.
No band-pass filter is involved, so a single volley counts as 1 cycle (a
10-15 Hz filter turns it into ~3 cycles of ripple, Edu-questions.md Q5).
`--self-test` checks the measure on a synthetic spike train.

Usage:
    python3 neuron/volley_stats.py --trains "Option 2=res/2026-09-07/ctx_nrn_45453403.npz"
    python3 neuron/volley_stats.py --self-test

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


# -- Stage 0: spike-based RE volley trains ----------------------------------
TRAIN_MIN_IVI = 60.0    # ms; volleys closer than this are one volley (<= 16.7 Hz)
TRAIN_MAX_IVI = 150.0   # ms; volleys further apart start a new train (>= 6.7 Hz)
SPINDLE_MIN_CYCLES = 6  # "real spindle": >= 6 cycles (PLAN-spindels.md goal table)


def re_volleys(t, g, R, tstop, min_part=0.10, skip_ms=0.0, smooth_ms=3.0):
    """RE population volleys from spikes: peaks of the RE rate at least
    TRAIN_MIN_IVI apart in which >= min_part of RE cells fire within +-20 ms.

    Returns a list of dicts: t (peak, ms), re_part / tc_part (fraction of
    cells firing: RE within +-20 ms of the peak, TC within [-40, +10] ms, i.e.
    the TC activity that leads into this RE volley), re_sd (SD of the RE
    cells' first-spike times in the volley, ms)."""
    rlo, rhi = R["re"]
    tlo, thi = R["tc"]
    mre = (g >= rlo) & (g < rhi)
    tre, gre = t[mre], g[mre]
    mtc = (g >= tlo) & (g < thi)
    ttc, gtc = t[mtc], g[mtc]
    r, b = A._rate(tre, tstop, 1.0, smooth_ms)
    r[b < skip_ms] = 0.0
    pk, _ = find_peaks(r, height=1e-9, distance=int(TRAIN_MIN_IVI))
    out = []
    for i in pk:
        p = b[i]
        m = (tre >= p - 20.0) & (tre < p + 20.0)
        cells = np.unique(gre[m])
        part = len(cells) / (rhi - rlo)
        if part < min_part:
            continue
        first = [tre[m][gre[m] == c].min() for c in cells]
        mt = (ttc >= p - 40.0) & (ttc < p + 10.0)
        out.append({"t": float(p), "re_part": part,
                    "tc_part": len(np.unique(gtc[mt])) / (thi - tlo),
                    "re_sd": float(np.std(first))})
    return out


def volley_trains(volleys):
    """Chain volleys whose interval is <= TRAIN_MAX_IVI into trains."""
    trains, cur = [], []
    for v in volleys:
        if cur and v["t"] - cur[-1]["t"] > TRAIN_MAX_IVI:
            trains.append(cur)
            cur = []
        cur.append(v)
    if cur:
        trains.append(cur)
    out = []
    for tr in trains:
        ts = np.array([v["t"] for v in tr])
        dur = float(ts[-1] - ts[0])
        part = np.array([v["re_part"] for v in tr])
        out.append({"t0": float(ts[0]), "t1": float(ts[-1]), "cycles": len(tr),
                    "duration": dur,
                    "freq": (len(tr) - 1) * 1000.0 / dur if len(tr) > 1 else float("nan"),
                    "re_part": part, "tc_part": np.array([v["tc_part"] for v in tr]),
                    "re_sd": np.array([v["re_sd"] for v in tr]),
                    # waxing/waning: participation peaks strictly inside the train
                    "wax_wane": bool(len(tr) >= 3 and 0 < int(np.argmax(part)) < len(tr) - 1)})
    return out


def train_summary(trains, span_ms):
    """Per-run summary of volley trains over span_ms of analysed time."""
    if not trains:
        return {"n_trains": 0}
    cyc = np.array([t["cycles"] for t in trains])
    starts = np.array([t["t0"] for t in trains])
    sp = [t for t in trains if t["cycles"] >= SPINDLE_MIN_CYCLES]
    multi = [t for t in trains if t["cycles"] > 1]
    return {
        "n_trains": len(trains),
        "trains_per_min": len(trains) / (span_ms / 60000.0),
        "cycles_median": float(np.median(cyc)), "cycles_max": int(cyc.max()),
        "frac_single": float(np.mean(cyc == 1)),
        "cycles_hist": {k: int((cyc == k).sum()) for k in range(1, min(cyc.max(), 12) + 1)},
        "n_spindles": len(sp), "spindles_per_min": len(sp) / (span_ms / 60000.0),
        "freq_median": float(np.median([t["freq"] for t in multi])) if multi else float("nan"),
        "duration_max": float(max(t["duration"] for t in trains)),
        "iti_median": float(np.median(np.diff(starts))) if len(starts) > 1 else float("nan"),
        "tc_part_median": float(np.median(np.concatenate([t["tc_part"] for t in trains]))),
        "re_part_median": float(np.median(np.concatenate([t["re_part"] for t in trains]))),
        "re_sd_median": float(np.median(np.concatenate([t["re_sd"] for t in trains]))),
    }


def run_trains(npz, skip_ms=1000.0):
    t, g, R, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    trains = volley_trains(re_volleys(t, g, R, tstop, skip_ms=skip_ms))
    return trains, train_summary(trains, tstop - skip_ms)


def print_train_summary(label, s):
    print(f"\n== {label}")
    if not s["n_trains"]:
        print("  no RE volleys")
        return
    print(f"  RE volley trains: {s['n_trains']} ({s['trains_per_min']:.1f}/min), "
          f"median inter-train interval {s['iti_median']:.0f} ms")
    print(f"  cycles per train: median {s['cycles_median']:.0f}, max {s['cycles_max']}, "
          f"single-volley trains {100 * s['frac_single']:.0f}%; histogram "
          + " ".join(f"{k}:{v}" for k, v in s["cycles_hist"].items()))
    print(f"  spindle-like trains (>= {SPINDLE_MIN_CYCLES} cycles): {s['n_spindles']} "
          f"({s['spindles_per_min']:.1f}/min); multi-cycle frequency median "
          f"{s['freq_median']:.1f} Hz; longest train {s['duration_max']:.0f} ms")
    print(f"  per cycle: RE participation {100 * s['re_part_median']:.0f}%, "
          f"TC participation {100 * s['tc_part_median']:.0f}%, "
          f"RE first-spike SD {s['re_sd_median']:.2f} ms")


def self_test():
    """Synthetic run: an 8-cycle 12 Hz train with a waxing/waning envelope,
    a single volley, and a 3-cycle 10 Hz train; plus sparse background."""
    rng = np.random.default_rng(0)
    R = {"tc": (0, 100), "re": (100, 130)}
    ts, gs = [], []

    def volley(t0, part_re, part_tc):
        for c in rng.choice(np.arange(100, 130), int(30 * part_re), replace=False):
            for k in range(3):   # a 3-spike RE burst
                ts.append(t0 + rng.normal(0, 2) + 3 * k); gs.append(c)
        for c in rng.choice(np.arange(0, 100), int(100 * part_tc), replace=False):
            ts.append(t0 - 25 + rng.normal(0, 3)); gs.append(c)

    env = [0.3, 0.5, 0.8, 1.0, 1.0, 0.8, 0.5, 0.3]
    for k, a in enumerate(env):
        volley(2000 + k * 1000 / 12.0, a, 0.4 * a)
    volley(5000, 0.9, 0.5)
    for k in range(3):
        volley(8000 + k * 100.0, 0.7, 0.3)
    bg = rng.uniform(0, 10000, 40)
    ts += list(bg); gs += list(rng.integers(0, 130, 40))
    npz = {"times": np.array(ts), "gids": np.array(gs), "ranges": np.array(R, dtype=object),
           "tstop": 10000.0}
    trains, s = run_trains(npz, skip_ms=0.0)
    got = [(tr["cycles"], round(tr["freq"], 1) if tr["cycles"] > 1 else None, tr["wax_wane"])
           for tr in trains]
    want = [(8, 12.0, True), (1, None, False), (3, 10.0, False)]
    ok = got == want
    print("[self-test] trains found (cycles, Hz, waxing/waning):", got)
    print("[self-test] expected:                               ", want)
    print("[self-test]", "PASS" if ok else "FAIL")
    return ok


def _thal_lfp(npz, fs=1000.0):
    """Thalamic LFP proxy: mean per-cell synaptic-kernel signal of TC and RE
    (ctx_analyze._synaptic_lfp), NOT z-scored, so runs with different
    activity are on one scale. Used for the isolated thalamus (no cortex)."""
    t, g, R, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    out = 0.0
    for k in THAL:
        lo, hi = R[k]
        lfp, bins = A._synaptic_lfp(t[(g >= lo) & (g < hi)], tstop, fs=fs)
        out = out + lfp / max(1, hi - lo) / len(THAL)
    return out - out[bins >= min(1000.0, tstop / 4)].mean(), bins


def plot_train_reconstruction(cases, out_png, window, lfp="cortex", title=None):
    """Literature-style SO-band + 10-15 Hz reconstruction (ctx_analyze's
    format) per case, shared y-scales, with BOTH event definitions:
    grey shading = band-pass events (ctx_analyze._detect_spindles), markers
    at the top = spike-based RE volley trains (Stage 0): black tick = single
    volley, orange = 2-5 cycles, green = >= 6 cycles (spindle-like).
    Blue lines = stimulus times (npz 'stim_t' or 'kick_t').
    cases: list of (label, npz, colour)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fs = 1000.0
    comp = []
    for label, npz, color in cases:
        if lfp == "cortex":
            _, x, _, bins, _ = A._composite_lfp(npz, fs=fs)
        else:
            x, bins = _thal_lfp(npz, fs)
        so = A._bandpass(x, fs, 0.5, 2.0)
        sp = A._bandpass(x, fs, 10.0, 15.0)
        st, en = A._detect_spindles(sp, bins, skip_ms=min(1000.0, float(npz["tstop"]) / 4))
        trains = volley_trains(re_volleys(npz["times"], npz["gids"], npz["ranges"].item(),
                                          float(npz["tstop"]), skip_ms=0.0))
        stim = []
        if "stim_t" in npz:
            stim = list(np.asarray(npz["stim_t"]))
        elif "kick_t" in npz:
            stim = [float(npz["kick_t"])]
        comp.append((label, color, bins, so + sp, sp, st, en, trains, stim))
    w0 = lambda b: (b >= window[0]) & (b < window[1])
    rmax = max(np.abs(c[3][w0(c[2])]).max() for c in comp) or 1.0
    smax = max(np.abs(c[4][w0(c[2])]).max() for c in comp) or 1.0
    fig, axes = plt.subplots(2 * len(comp), 1, sharex=True,
                             figsize=(14, 2.25 * 2 * len(comp) + 0.9))
    for i, (label, color, bins, rec, sp, st, en, trains, stim) in enumerate(comp):
        w = w0(bins)
        ts = bins[w] / 1000.0
        a1, a2 = axes[2 * i], axes[2 * i + 1]
        a1.plot(ts, rec[w], color=color, lw=0.8)
        a1.set_ylim(-rmax * 1.1, rmax * 1.35)
        a1.set_ylabel("raw\n(SO+spindle)", fontsize=8)
        a1.set_title(label, fontsize=10, loc="left", fontweight="bold")
        a2.plot(ts, sp[w], color=color, lw=0.8)
        a2.set_ylim(-smax * 1.1, smax * 1.1)
        a2.set_ylabel("10-15 Hz", fontsize=8)
        for s_, e_ in zip(st, en):
            a, b = bins[s_] / 1000.0, bins[e_] / 1000.0
            if b >= window[0] / 1000.0 and a <= window[1] / 1000.0:
                for ax in (a1, a2):
                    ax.axvspan(a, b, color="0.5", alpha=0.18, lw=0)
        y = rmax * 1.22
        for tr in trains:
            if not (window[0] <= tr["t0"] < window[1]):
                continue
            c = tr["cycles"]
            col = "k" if c == 1 else ("#e67e22" if c < SPINDLE_MIN_CYCLES else "#1e8449")
            if c == 1:
                a1.plot([tr["t0"] / 1000.0] * 2, [y - rmax * 0.07, y + rmax * 0.07],
                        color=col, lw=1.2)
            else:
                a1.plot([tr["t0"] / 1000.0, tr["t1"] / 1000.0], [y, y], color=col, lw=4,
                        solid_capstyle="butt")
                a1.text(tr["t1"] / 1000.0 + 0.02, y, str(c), color=col, fontsize=7,
                        va="center")
        for t in stim:
            if window[0] <= t < window[1]:
                for ax in (a1, a2):
                    ax.axvline(t / 1000.0, color=A.BLUE, lw=1.2, alpha=0.8)
    axes[-1].set_xlabel("time (s)")
    axes[-1].set_xlim(window[0] / 1000.0, window[1] / 1000.0)
    src = ("cortical LFP proxy (L2/3-L6 composite)" if lfp == "cortex"
           else "thalamic LFP proxy (TC+RE, per cell)")
    import textwrap
    cap = textwrap.fill(src + ", SO band (0.5-2 Hz) + spindle band (10-15 Hz). Grey = "
                        "band-pass events; top markers = spike-based RE volley trains "
                        "(black tick = 1 cycle, orange = 2-5, green = >= 6 cycles, "
                        "number = cycles); blue = stimulus", 150)
    fig.suptitle((title + "\n" if title else "") + cap, fontsize=9)
    top = 1 - (0.05 + 0.17 * (cap.count("\n") + (2 if title else 1))) / fig.get_size_inches()[1]
    fig.tight_layout(rect=[0, 0, 1, top])
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    print(f"Saved {out_png}")


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
    ap.add_argument("runs", nargs="*", help='"label=path.npz" pairs')
    ap.add_argument("--trains", action="store_true",
                    help="Stage 0 spike-based RE volley-train report instead")
    ap.add_argument("--recon", default="",
                    help="output PNG: literature-style SO+10-15 Hz reconstruction "
                         "per run with band-pass events AND spike-based trains "
                         "(uses --window-start/--window-len, --lfp)")
    ap.add_argument("--lfp", choices=["cortex", "thal"], default="cortex")
    ap.add_argument("--title", default=None)
    ap.add_argument("--self-test", action="store_true",
                    help="check the train measure on a synthetic spike train")
    ap.add_argument("--plot", default=None, help="optional output PNG for the rate plot")
    ap.add_argument("--window-start", type=float, default=20000.0)
    ap.add_argument("--window-len", type=float, default=3000.0)
    ap.add_argument("--windows", default="",
                    help='per-window table instead of the full-run report, e.g. '
                         '"20-40,60-80,100-120,160-180" (seconds)')
    a = ap.parse_args(argv)

    if a.self_test:
        return 0 if self_test() else 1
    cases = []
    for pair in a.runs:
        label, path = pair.split("=", 1)
        cases.append((label, np.load(path, allow_pickle=True)))
    if a.recon:
        cols = ["#b03a2e", "0.3", "#1f618d", "#7d3c98", "#117a65", "#b9770e"]
        plot_train_reconstruction([(l, n, cols[i % len(cols)]) for i, (l, n) in enumerate(cases)],
                                  a.recon, (a.window_start, a.window_start + a.window_len),
                                  lfp=a.lfp, title=a.title)
        return 0
    if a.trains:
        for label, npz in cases:
            print_train_summary(label, run_trains(npz)[1])
        return 0
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
