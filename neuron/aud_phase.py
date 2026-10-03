"""
Auditory plan Stages C1 / C2 (PLAN-auditory-input.md): how a tone's effect
depends on the slow-oscillation (SO) phase and on the spindle state at the
moment it arrives. No NEURON needed.

Each tone run is paired with a sham run: the same network (same flags and
HET_SEED, so the same background-input trains) without --stim. The tone
times of the stim run are applied to the sham run as "fake tones", and both
runs are classified and measured the same way, each on its own activity. The
difference stim - sham per bin is the tone's effect in that state.

State at tone onset t0, from activity BEFORE t0 only:
  - SO phase (C1), from the L5E + L6E rate in [t0 - 100, t0):
      UP          rate >= --up-thr (Hz/cell; the rate is bimodal, DOWN ~0)
      DOWN-early  below threshold, time since the last UP onset < the median
                  of that time over DOWN tones (pooled over runs)
      DOWN-late   below threshold, later (closer to the next spontaneous UP)
  - spindle state (C2), from RE volley trains (volley_stats.py) that started
    before t0:
      in-spindle    t0 within a >= 3-cycle train (first .. last volley + 50 ms)
      post-spindle  0-2 s after the end of a >= 6-cycle train (refractory)
      none          otherwise

Responses per tone:
  - TC, L4E, L2/3E, L5E, L6E spikes/cell in 0-50 ms; TC fraction of those
    spikes in Lu-criterion bursts (aud_analyze.burst_mask)
  - evoked UP: an UP onset (stage2_report.up_onsets) in (t0, t0 + 300 ms]
  - evoked spindle: a >= 6-cycle RE train starting in (t0, t0 + 300 ms]

    python3 neuron/aud_phase.py --pair stim.npz=sham.npz [--pair ...] \\
        [--out-txt table.txt] [--plot fig.png]

Several --pair arguments (e.g. seeds) are pooled.

Without a sham run (--pair stim.npz, no "="), the control is provisional
and comes from the same run: one control time per inter-tone interval, drawn
uniformly in [previous tone + 1.5 s, next tone - 0.5 s] (deterministic). Those
times are still within a few seconds of a tone, so slow after-effects of the
previous tone (e.g. spindle refractoriness) are not removed; use a sham run
for the real comparison.
"""

import argparse
import os
import sys

import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aud_analyze as AA                                 # noqa: E402
import ctx_analyze as A                                  # noqa: E402
import stage2_report as S                                # noqa: E402
import volley_stats as V                                 # noqa: E402

POPS = ("tc", "l4e", "l23e", "l5e", "l6e")
SO_BINS = ("UP", "DOWN-early", "DOWN-late")
SP_BINS = ("in-spindle", "post-spindle", "none")


def measure(npz, tones, skip=2000.0, up_thr=2.0):
    """Per-tone state and responses for one run. tones: onset times (ms)."""
    t, g, R, T = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    tones = np.asarray(tones, float)
    tones = tones[(tones >= skip) & (tones < T - 300.0)]
    ups = S.up_onsets(t, g, R, T, skip)
    trains = V.volley_trains(V.re_volleys(t, g, R, T, skip_ms=skip))

    m = np.zeros_like(t, bool)
    n_cx = 0
    for k in ("l5e", "l6e"):
        lo, hi = R[k]
        m |= (g >= lo) & (g < hi)
        n_cx += hi - lo
    tcx = np.sort(t[m])

    per_pop = {}
    for p in POPS:
        lo, hi = R[p]
        sel = (g >= lo) & (g < hi)
        per_pop[p] = (np.sort(t[sel]), hi - lo)
    lo, hi = R["tc"]
    sel = (g >= lo) & (g < hi)
    ttc = t[sel]
    burst = AA.burst_mask(ttc, g[sel])
    ttc_b = np.sort(ttc[burst])

    def count(ts, a, b):
        return np.searchsorted(ts, b) - np.searchsorted(ts, a)

    rows = []
    for t0 in tones:
        pre = count(tcx, t0 - 100.0, t0) / n_cx / 0.1
        prev = ups[ups < t0]
        since = t0 - prev[-1] if len(prev) else np.nan
        so = "UP" if pre >= up_thr else "DOWN"
        sp = "none"
        for tr in trains:
            if tr["t0"] >= t0:
                break
            if tr["cycles"] >= 3 and tr["t0"] <= t0 <= tr["t1"] + 50.0:
                sp = "in-spindle"
            elif tr["cycles"] >= 6 and 0 < t0 - tr["t1"] <= 2000.0 and sp == "none":
                sp = "post-spindle"
        row = {"t0": t0, "pre_rate": pre, "since_up": since, "so": so, "sp": sp}
        for p in POPS:
            ts, n = per_pop[p]
            row[p] = count(ts, t0, t0 + 50.0) / n
        n_tc_sp = count(per_pop["tc"][0], t0, t0 + 50.0)
        row["tc_burst"] = count(ttc_b, t0, t0 + 50.0) / n_tc_sp if n_tc_sp else np.nan
        row["evoked_up"] = bool(np.any((ups > t0) & (ups <= t0 + 300.0)))
        row["evoked_sp"] = any(t0 < tr["t0"] <= t0 + 300.0 and tr["cycles"] >= 6
                               for tr in trains)
        rows.append(row)
    return rows


def control_times(tones, seed=0):
    """Within-run control times between consecutive tones (see docstring)."""
    rng = np.random.default_rng(seed)
    out, src = [], []
    for i in range(len(tones) - 1):
        lo, hi = tones[i] + 1500.0, tones[i + 1] - 500.0
        if hi > lo:
            out.append(rng.uniform(lo, hi))
            src.append(i)
    return np.array(out), np.array(src, int)


def split_down(rows_all):
    """DOWN -> DOWN-early / DOWN-late at the pooled median time since UP."""
    since = [r["since_up"] for r in rows_all if r["so"] == "DOWN" and np.isfinite(r["since_up"])]
    cut = float(np.median(since)) if since else np.inf
    for r in rows_all:
        if r["so"] == "DOWN":
            r["so"] = "DOWN-early" if (np.isfinite(r["since_up"]) and r["since_up"] < cut) else "DOWN-late"
    return cut


def _cell(stim, sham, key):
    a = np.array([r[key] for r in stim], float)
    b = np.array([r[key] for r in sham], float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if not len(a) or not len(b):
        return "%-24s" % "-", np.nan, np.nan
    if key in ("evoked_up", "evoked_sp"):
        p = fisher_exact([[a.sum(), len(a) - a.sum()], [b.sum(), len(b) - b.sum()]])[1]
        return "%3.0f%% vs %3.0f%% (p %.2g)" % (100 * a.mean(), 100 * b.mean(), p), a.mean(), b.mean()
    p = mannwhitneyu(a, b).pvalue if len(a) > 1 and len(b) > 1 else np.nan
    return "%5.2f vs %5.2f (p %.2g)" % (a.mean(), b.mean(), p), a.mean(), b.mean()


KEYS = [("tc", "TC sp/cell 0-50"), ("tc_burst", "TC burst frac"),
        ("l4e", "L4E sp/cell"), ("l23e", "L2/3E sp/cell"),
        ("l5e", "L5E sp/cell"), ("l6e", "L6E sp/cell"),
        ("evoked_up", "evoked UP (300 ms)"), ("evoked_sp", "evoked spindle")]


def table(stim, sham, field, bins, title, levels=None):
    lines = ["", "== %s: tone vs sham, mean (Mann-Whitney / Fisher p)" % title]
    for lev in (levels if levels is not None else [None]):
        st = [r for r in stim if lev is None or r["level"] == lev]
        sh = [r for r in sham if lev is None or r["level"] == lev]
        if lev is not None:
            lines.append("  -- level %g dB" % lev)
        hdr = "  %-20s" % "" + "".join("| %-26s" % b for b in bins)
        lines.append(hdr)
        lines.append("  %-20s" % "n tones (tone/sham)" + "".join(
            "| %-26s" % ("%d / %d" % (sum(r[field] == b for r in st), sum(r[field] == b for r in sh)))
            for b in bins))
        for key, name in KEYS:
            cells = [_cell([r for r in st if r[field] == b], [r for r in sh if r[field] == b], key)[0]
                     for b in bins]
            lines.append("  %-20s" % name + "".join("| %-26s" % c for c in cells))
    return lines


def figure(stim, sham, out_png, ctl_label="sham"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [("so", SO_BINS, "SO phase at tone onset (C1)"),
              ("sp", SP_BINS, "spindle state at tone onset (C2)")]
    keys = [("tc", "TC spikes/cell, 0-50 ms"), ("l4e", "L4E spikes/cell, 0-50 ms"),
            ("evoked_up", "P(UP onset within 300 ms)"), ("evoked_sp", "P(spindle within 300 ms)")]
    fig, axes = plt.subplots(len(panels), len(keys), figsize=(15, 6.5))
    for i, (field, bins, title) in enumerate(panels):
        for j, (key, name) in enumerate(keys):
            ax = axes[i, j]
            x = np.arange(len(bins))
            vs = [_cell([r for r in stim if r[field] == b], [r for r in sham if r[field] == b], key)
                  for b in bins]
            ax.bar(x - 0.2, [v[1] for v in vs], 0.4, color="#b03a2e", label="tone")
            ax.bar(x + 0.2, [v[2] for v in vs], 0.4, color="0.6", label=ctl_label)
            ax.set_xticks(x)
            ax.set_xticklabels(["%s\n(n=%d)" % (b, sum(r[field] == b for r in stim)) for b in bins],
                               fontsize=8)
            ax.set_title(name, fontsize=9)
            if j == 0:
                ax.set_ylabel(title, fontsize=9)
            if i == 0 and j == 0:
                ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    print("Saved", out_png)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", action="append", required=True,
                    help="stim.npz=sham.npz (repeat to pool seeds)")
    ap.add_argument("--skip", type=float, default=2000.0)
    ap.add_argument("--up-thr", type=float, default=2.0,
                    help="UP if the L5E+L6E rate in the 100 ms before the tone >= this (Hz/cell)")
    ap.add_argument("--out-txt", default="")
    ap.add_argument("--plot", default="")
    a = ap.parse_args(argv)

    stim_all, sham_all = [], []
    self_ctl = False
    for pair in a.pair:
        sp, sh = pair.split("=", 1) if "=" in pair else (pair, None)
        zs = np.load(sp, allow_pickle=True)
        if "stim_t" not in zs:
            raise SystemExit("%s has no stimulus log (run without --stim?)" % sp)
        tones, levels = zs["stim_t"], zs["stim_level"]
        lev_of = dict(zip(np.round(tones, 6), levels))
        for r in measure(zs, tones, a.skip, a.up_thr):
            r["level"] = float(lev_of[round(r["t0"], 6)])
            r["run"] = os.path.basename(sp)
            stim_all.append(r)
        if sh is not None:
            ctl, ctl_lev = tones, lev_of
            z = np.load(sh, allow_pickle=True)
        else:
            self_ctl = True
            ctl, src = control_times(np.asarray(tones, float))
            ctl_lev = dict(zip(np.round(ctl, 6), np.asarray(levels)[src]))
            z = zs
        for r in measure(z, ctl, a.skip, a.up_thr):
            r["level"] = float(ctl_lev[round(r["t0"], 6)])
            r["run"] = os.path.basename(sh or sp) + (" (within-run control)" if sh is None else "")
            sham_all.append(r)
    cut = split_down(stim_all + sham_all)
    levels = sorted({r["level"] for r in stim_all})
    out = (["PROVISIONAL: control = within-run times between tones (no sham run); "
            "'sham' columns below are those controls"] if self_ctl else []) + ["Auditory C1/C2: %d tone runs, %d tones (levels %s dB); DOWN-early/late split at "
           "%.0f ms after the last UP onset; UP threshold %.1f Hz/cell"
           % (len(a.pair), len(stim_all), ", ".join("%g" % x for x in levels), cut, a.up_thr)]
    out += table(stim_all, sham_all, "so", SO_BINS, "C1, SO phase (all levels)")
    out += table(stim_all, sham_all, "sp", SP_BINS, "C2, spindle state (all levels)")
    if len(levels) > 1:
        out += table(stim_all, sham_all, "so", SO_BINS, "C1, SO phase, by level", levels)
    text = "\n".join(out)
    print(text)
    if a.out_txt:
        with open(a.out_txt, "w") as f:
            f.write(text + "\n")
    if a.plot:
        figure(stim_all, sham_all, a.plot,
               "within-run control (provisional)" if self_ctl else "sham")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
