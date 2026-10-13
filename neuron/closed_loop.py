"""
C3 closed-loop stimulation, step 1: a causal online slow-oscillation (SO)
state detector, and an offline replay that calibrates it on recorded spikes.

PLAN-auditory-input.md D5: run psolve in chunks, allreduce the L5E + L6E
spike count after each chunk, feed a causal estimator, stimulate at a delay d
after a DOWN state is detected (Ngo et al. 2013), with a mandatory sham
(same detector, detections logged, no sound).

The detector is deliberately a plain class that sees only (chunk end time,
spike count). The online driver and the offline replay below feed it the same
way, so what is calibrated offline is exactly what runs in the simulation.
Every rank gets the same allreduced counts and runs the same deterministic
detector, so no decision needs broadcasting.

State definitions follow aud_phase.py: UP = L5E+L6E rate over the last 100 ms
>= up_thr (2.0 Hz/cell; the rate is bimodal, DOWN ~ 0).

    python3 neuron/closed_loop.py res/2026-10-09/c1c2/c12_sham_s0.npz
"""

import argparse
from collections import deque

import numpy as np


class OnlineSODetector:
    """Causal detector of the slow-oscillation UP volley on a stream of
    per-chunk L5E+L6E spike counts.

    Why UP-onset and not DOWN: in this model the cortex never goes silent.
    Between slow oscillations the L5E+L6E rate is 2-3.5 Hz/cell (background
    input), troughs of the 100 ms rate reach only ~0.9, while the UP onset is
    a sharp ~10 Hz/cell volley. A "population silence" DOWN detector (the
    PLAN's wording) has nothing to detect; the volley is the robust landmark.
    Ngo et al. (2013) also stimulate at a fixed delay after the detected
    up-state, so the phase logic is the same.

    step(t_end, count) returns [("up", t_end)] when the windowed rate first
    crosses up_thr while armed. t_end is the END of the chunk, i.e. the
    earliest time a causal observer could know. After a detection the
    detector re-arms once refractory_ms has passed and the rate is below
    up_thr / 2.
    """

    def __init__(self, n_cx, chunk_ms=5.0, win_ms=50.0, up_thr=6.0, refractory_ms=300.0):
        self.n_cx = int(n_cx)
        self.chunk_ms = float(chunk_ms)
        self.n_win = max(1, int(round(win_ms / chunk_ms)))
        self.win_s = self.n_win * chunk_ms / 1000.0
        self.up_thr = float(up_thr)
        self.refractory_ms = float(refractory_ms)
        self.buf = deque([0] * self.n_win, maxlen=self.n_win)
        self.armed = True
        self.t_last = -1e18

    def rate(self):
        """Window rate, Hz per cell."""
        return sum(self.buf) / self.n_cx / self.win_s

    def step(self, t_end, count):
        self.buf.append(int(count))
        r = self.rate()
        if (not self.armed and t_end - self.t_last >= self.refractory_ms
                and r < 0.5 * self.up_thr):
            self.armed = True
        if self.armed and r >= self.up_thr:
            self.armed = False
            self.t_last = float(t_end)
            return [("up", float(t_end))]
        return []


def run_closed_loop(net, aud, mode, delays, every=2, level=60.0, dur=50.0,
                    seed=1, chunk_ms=5.0, tstop=None, celsius=36.0, dt=0.025,
                    win_ms=50.0, up_thr=6.0, skip=2000.0):
    """Drop-in replacement for net.run() that runs psolve in chunks and, on
    each detected UP volley, queues a tone at a delay d (mode "stim") or does
    nothing (mode "sham": same detector, same delay draws, detections logged,
    no sound -- the matched control the plan requires).

    d is drawn per detection from `delays` with an RNG seeded by (seed,
    detection index), so a stim run and its sham make identical draws and
    can be compared at matched detection times. Only every `every`-th
    detection (after `skip`) is stimulated, to keep the tone spacing
    comparable to the open-loop runs (2-5 s).

    Every rank sees the same allreduced count, runs the same deterministic
    detector, and so makes the same decision without any broadcast.
    Returns (wall_seconds, log_dict).
    """
    import time
    from neuron import h
    assert mode in ("stim", "sham")
    if mode == "stim":
        assert aud is not None, "stim mode needs an AuditoryInput"
    pc = net.pc
    pc.set_maxstep(10)
    h.celsius, h.dt = celsius, dt
    h.finitialize(-70)

    cx = [net.ranges[k] for k in ("l5e", "l6e")]
    n_cx = sum(hi - lo for lo, hi in cx)
    det = OnlineSODetector(n_cx, chunk_ms, win_ms, up_thr)
    log = {"t_det": [], "d": [], "stim": [], "level": []}
    ptr, k_det, n_stimmed = 0, 0, 0
    t0 = time.time()
    n_chunks = int(round(tstop / chunk_ms))
    for i in range(n_chunks):
        t_end = (i + 1) * chunk_ms
        pc.psolve(t_end)
        gs = np.asarray(net.gspk)
        new = gs[ptr:]
        ptr = len(gs)
        local = int(sum(((new >= lo) & (new < hi)).sum() for lo, hi in cx))
        total = int(pc.allreduce(local, 1))
        for _, t_det in det.step(t_end, total):
            if t_det < skip:
                continue
            r = np.random.default_rng([int(seed), 11, k_det])
            d = float(r.choice(delays))
            stimmed = (k_det % every == 0)
            lv = float(np.atleast_1d(level)[r.integers(len(np.atleast_1d(level)))])
            log["t_det"].append(t_det); log["d"].append(d)
            log["stim"].append(int(stimmed)); log["level"].append(lv)
            if mode == "stim" and stimmed and t_det + d < tstop - 400.0:
                aud.queue_tone(t_det + d, dur, lv)
                n_stimmed += 1
            k_det += 1
    wall = pc.allreduce(time.time() - t0, 2)
    out = {"cl_mode": mode, "cl_t_det": np.array(log["t_det"]), "cl_d": np.array(log["d"]),
           "cl_stim": np.array(log["stim"], int), "cl_level": np.array(log["level"]),
           "cl_params": "win_ms=%g up_thr=%g chunk_ms=%g every=%d delays=%s"
                        % (win_ms, up_thr, chunk_ms, every, list(delays))}
    return wall, out


def chunk_counts(t, g, ranges, tstop, chunk_ms=5.0):
    """Per-chunk L5E+L6E spike counts from recorded spikes, and n_cx."""
    m = np.zeros(len(t), bool)
    n_cx = 0
    for k in ("l5e", "l6e"):
        lo, hi = ranges[k]
        m |= (g >= lo) & (g < hi)
        n_cx += hi - lo
    edges = np.arange(0.0, tstop + chunk_ms, chunk_ms)
    return np.histogram(t[m], bins=edges)[0], n_cx


def replay(npz, chunk_ms=5.0, win_ms=50.0, up_thr=6.0, refractory_ms=300.0, skip=2000.0):
    """Run the detector over a recorded run. Returns detection times (ms)."""
    t, g = npz["times"], npz["gids"]
    R, T = npz["ranges"].item(), float(npz["tstop"])
    counts, n_cx = chunk_counts(t, g, R, T, chunk_ms)
    det = OnlineSODetector(n_cx, chunk_ms, win_ms, up_thr, refractory_ms)
    out = []
    for i, c in enumerate(counts):
        for _, tt in det.step((i + 1) * chunk_ms, c):
            if tt >= skip:
                out.append(tt)
    return np.array(out)


def score(det_t, ups, tol=150.0):
    """Match detections to offline UP onsets (stage2_report.up_onsets, the
    peak of the 20 ms rate). A detection within [-tol, +tol] of an onset is a
    hit. Returns hit rate, false-positive fraction, latency (det - onset)."""
    ups = np.asarray(ups)
    lat, used = [], set()
    for td in det_t:
        j = int(np.argmin(np.abs(ups - td))) if len(ups) else -1
        if j >= 0 and abs(ups[j] - td) <= tol and j not in used:
            used.add(j)
            lat.append(td - ups[j])
    n_hit = len(lat)
    return (n_hit / max(len(ups), 1), 1 - n_hit / max(len(det_t), 1), np.array(lat))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("npz", nargs="+")
    ap.add_argument("--win-ms", type=float, default=50.0)
    ap.add_argument("--up-thr", type=float, default=6.0)
    ap.add_argument("--sweep", action="store_true",
                    help="sweep window and threshold instead of the phase table")
    ap.add_argument("--delays", default="0,200,400,600,800,1000,1200,1400")
    ap.add_argument("--down-cut", type=float, default=1273.0,
                    help="DOWN-early/late split, ms since the last UP onset "
                         "(C1/C2 final table)")
    args = ap.parse_args(argv)

    import stage2_report as S    # offline ground truth

    runs = []
    for path in args.npz:
        npz = np.load(path, allow_pickle=True)
        t, g = npz["times"], npz["gids"]
        R, T = npz["ranges"].item(), float(npz["tstop"])
        runs.append((npz, S.up_onsets(t, g, R, T, 2000.0), T))
    n_ups = sum(len(u) for _, u, _ in runs)

    if args.sweep:
        print(f"{len(runs)} runs, {n_ups} offline UP onsets")
        print(f"{'win ms':>7} {'thr Hz/cell':>12} | {'hit rate':>9} {'false-pos':>10} | "
              f"{'latency median (IQR) ms':>26}")
        for win in (10, 20, 50):
            for thr in (4, 5, 6, 8, 10):
                hits, fps, lats, nd = [], [], [], 0
                for npz, ups, T in runs:
                    d = replay(npz, win_ms=win, up_thr=thr)
                    h, f, l = score(d, ups)
                    nd += len(d); hits.append(h * len(ups)); lats.append(l)
                    fps.append(f * len(d))
                l = np.concatenate(lats)
                print(f"{win:7d} {thr:12d} | {sum(hits)/n_ups:9.0%} {sum(fps)/max(nd,1):10.0%} | "
                      f"{np.median(l):10.0f} ({np.percentile(l,25):.0f}..{np.percentile(l,75):.0f})")
        return 0

    delays = [float(x) for x in args.delays.split(",")]
    rows = {d: [] for d in delays}
    lat_all, nd_all = [], 0
    for npz, ups, T in runs:
        t, g, R = npz["times"], npz["gids"], npz["ranges"].item()
        dets = replay(npz, win_ms=args.win_ms, up_thr=args.up_thr)
        nd_all += len(dets)
        _, _, l = score(dets, ups)
        lat_all.append(l)
        counts, n_cx = chunk_counts(t, g, R, T, 5.0)
        csum = np.concatenate([[0], np.cumsum(counts)])
        for td in dets:
            for d in delays:
                ts = td + d
                if ts >= T - 300.0:
                    continue
                i1 = int(ts / 5.0)
                pre = (csum[i1] - csum[max(0, i1 - 20)]) / n_cx / 0.1
                prev, nxt = ups[ups < ts], ups[ups > ts]
                since = ts - prev[-1] if len(prev) else np.nan
                state = "UP" if pre >= 2.0 else ("DOWN-early" if since < args.down_cut else "DOWN-late")
                rows[d].append((state, since, (nxt[0] - ts) if len(nxt) else np.nan))
    l = np.concatenate(lat_all)
    print(f"{len(runs)} runs; detector win {args.win_ms:.0f} ms, thr {args.up_thr:.1f} Hz/cell: "
          f"{nd_all} detections vs {n_ups} offline UP onsets")
    print(f"latency of detection vs offline UP onset: median {np.median(l):.0f} ms "
          f"(IQR {np.percentile(l,25):.0f}..{np.percentile(l,75):.0f})\n")
    print("tone placed at (detection + d): SO state at that moment (aud_phase definitions), "
          "and time to the next UP onset")
    print(f"{'d (ms)':>7} {'n':>5} | {'UP':>6} {'DOWN-early':>11} {'DOWN-late':>10} | "
          f"{'since UP, median (IQR) ms':>27} | {'to next UP, median (IQR) ms':>28}")
    for d in delays:
        r = rows[d]
        n = len(r)
        c = {k: sum(1 for x in r if x[0] == k) / n for k in ("UP", "DOWN-early", "DOWN-late")}
        su = np.array([x[1] for x in r if np.isfinite(x[1])])
        nu = np.array([x[2] for x in r if np.isfinite(x[2])])
        print(f"{d:7.0f} {n:5d} | {c['UP']:6.0%} {c['DOWN-early']:11.0%} {c['DOWN-late']:10.0%} | "
              f"{np.median(su):8.0f} ({np.percentile(su,25):.0f}..{np.percentile(su,75):.0f}) | "
              f"{np.median(nu):9.0f} ({np.percentile(nu,25):.0f}..{np.percentile(nu,75):.0f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
