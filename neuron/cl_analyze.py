"""
C3 closed-loop analysis: stim run vs its matched sham, binned by the delay d
after the detected UP volley (PLAN-auditory-input.md D5/C3).

Both arms of a pair run the same detector and make the same delay draws
(closed_loop.run_closed_loop); only the stim arm plays the tone. A stim run's
dynamics diverge from its sham after the first tone, so the arms are NOT
compared at identical detection times. They are compared the way C1/C2 did:
for every scheduled tone, measure the response after (detection + d), in the
stim run for real and in the sham run for the fake tone at the same d.

Per d (and pooled):
  P(UP onset within 300 ms), P(spindle within 300 ms), TC spikes/cell in 50 ms
  -- Fisher / Mann-Whitney p, stim vs sham.
Per arm, over the whole run (what the loop does to the slow oscillation):
  UP-onset rate, CV of UP intervals, spindles per UP, spindle rate.

    python3 neuron/cl_analyze.py --pair out/cl_stim_s0.npz=out/cl_sham_s0.npz \\
        --pair out/cl_stim_s1.npz=out/cl_sham_s1.npz --out-txt res/cl_table.txt
"""

import argparse

import numpy as np
from scipy import stats

import aud_phase as AP
import stage2_report as S
import volley_stats as V


def tones_of(npz, tstop_guard=400.0):
    """(tone times, d) of a closed-loop run's scheduled tones. For a sham run
    these are the fake tones; for a stim run they are the real ones. Same
    rule as run_closed_loop: every `cl_stim`-flagged detection, unless the
    tone would fall within tstop_guard of the end."""
    T = float(npz["tstop"])
    sel = npz["cl_stim"] == 1
    t = npz["cl_t_det"][sel] + npz["cl_d"][sel]
    d = npz["cl_d"][sel]
    ok = t < T - tstop_guard
    return t[ok], d[ok]


def so_stats(npz, skip=2000.0):
    t, g, R, T = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    ups = S.up_onsets(t, g, R, T, skip)
    trains = [tr for tr in V.volley_trains(V.re_volleys(t, g, R, T, skip_ms=skip))
              if tr["cycles"] >= 6]
    sp_t0 = np.array([tr["t0"] for tr in trains])
    iv = np.diff(ups)
    with_sp = sum(1 for u in ups if np.any((sp_t0 > u) & (sp_t0 <= u + 300.0)))
    dur = (T - skip) / 1000.0
    return {"up_hz": len(ups) / dur,
            "up_cv": float(iv.std() / iv.mean()) if len(iv) > 2 else np.nan,
            "sp_per_up": with_sp / max(len(ups), 1),
            "sp_hz": len(trains) / dur,
            "n_up": len(ups), "n_sp": len(trains)}


def collect(pairs, skip, up_thr):
    stim_rows, sham_rows, stim_so, sham_so = [], [], [], []
    for sp, sh in pairs:
        zs, zh = np.load(sp, allow_pickle=True), np.load(sh, allow_pickle=True)
        for z, rows, so in ((zs, stim_rows, stim_so), (zh, sham_rows, sham_so)):
            t, d = tones_of(z)
            for r, dd in zip(AP.measure(z, t, skip, up_thr), d):
                r["d"] = float(dd)
                rows.append(r)
            so.append(so_stats(z, skip))
    return stim_rows, sham_rows, stim_so, sham_so


def _p_bool(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    return stats.fisher_exact([[a.sum(), (~a).sum()], [b.sum(), (~b).sum()]])[1]


def _p_cont(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return stats.mannwhitneyu(a, b)[1]


def table(stim, sham, delays):
    lines = [f"{'d (ms)':>7} | {'n tone/sham':>11} | {'P(UP<300ms)':>22} | {'P(spindle<300ms)':>22} | "
             f"{'TC sp/cell (50ms)':>26}", "-" * 104]
    for d in list(delays) + ["all"]:
        a = [r for r in stim if d == "all" or r["d"] == d]
        b = [r for r in sham if d == "all" or r["d"] == d]
        if not a or not b:
            continue
        f = lambda rows, k: np.array([r[k] for r in rows])
        up_a, up_b = f(a, "evoked_up").mean(), f(b, "evoked_up").mean()
        sp_a, sp_b = f(a, "evoked_sp").mean(), f(b, "evoked_sp").mean()
        tc_a, tc_b = f(a, "tc").mean(), f(b, "tc").mean()
        lines.append(
            f"{str(d):>7} | {len(a):>4}/{len(b):<6} | "
            f"{up_a:5.0%} vs {up_b:4.0%} (p={_p_bool(f(a,'evoked_up'), f(b,'evoked_up')):.3f}) | "
            f"{sp_a:5.0%} vs {sp_b:4.0%} (p={_p_bool(f(a,'evoked_sp'), f(b,'evoked_sp')):.3f}) | "
            f"{tc_a:5.2f} vs {tc_b:4.2f} (p={_p_cont(f(a,'tc'), f(b,'tc')):.3f})")
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", action="append", required=True,
                    help="stim.npz=sham.npz (repeat per seed)")
    ap.add_argument("--skip", type=float, default=2000.0)
    ap.add_argument("--up-thr", type=float, default=2.0)
    ap.add_argument("--out-txt", default="")
    a = ap.parse_args(argv)

    pairs = [p.split("=", 1) for p in a.pair]
    stim, sham, so_stim, so_sham = collect(pairs, a.skip, a.up_thr)
    delays = sorted({r["d"] for r in stim + sham})

    out = [f"C3 closed loop: {len(pairs)} stim/sham pair(s); {len(stim)} tones (stim), "
           f"{len(sham)} fake tones (sham); response windows 300 ms (UP, spindle), 50 ms (TC)", ""]
    out += table(stim, sham, delays)
    out += ["", "Whole-run slow oscillation, per seed (stim | sham):",
            f"{'seed':>5} | {'UP Hz':>12} | {'UP CV':>12} | {'spindles/UP':>14} | {'spindle Hz':>14}"]
    for i, (x, y) in enumerate(zip(so_stim, so_sham)):
        out.append(f"{i:>5} | {x['up_hz']:5.2f} | {y['up_hz']:4.2f} | {x['up_cv']:5.2f} | {y['up_cv']:4.2f} | "
                   f"{x['sp_per_up']:6.2f} | {y['sp_per_up']:5.2f} | {x['sp_hz']:6.3f} | {y['sp_hz']:5.3f}")
    txt = "\n".join(out)
    print(txt)
    if a.out_txt:
        with open(a.out_txt, "w") as f:
            f.write(txt + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
