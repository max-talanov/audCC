"""
Can the isolated thalamus ring? (PLAN-spindels.md Stage 1)

Builds the thalamus of ParallelCorticoThalamicNet at MN5 size (346 TC / 91
RE, --scale 1.65) with every cortical population set to 0 cells. So the TC
and RE cells, their heterogeneity, RE->TC / TC->RE / RE<->RE wiring and gap
junctions are exactly the production ones, and the L6 -> TC / RE synapses
exist but have no presynaptic cells. (Not tc_mpi.py: it still defaults to the
old gsk_re = 5e-5 and uses a periodic drive.)

Protocol per run: 1 s settle, then ONE kick at --kick-t, then free-running to
--tstop (>= 3 s after the kick). Kick types:
  l6    a single L6 corticothalamic volley: each TC / RE cell gets 100 AMPA
        events on its production L6 synapse, total conductance g_l6_tc /
        g_l6_re (0.03 uS each, as in production), times jittered with SD
        --kick-sd ms around the kick time
  tone  one 60 dB, 50 ms tone through the IC -> TC driver input
        (auditory_input.py, Stage A calibration); RE is reached only via TC

Seeds vary the per-cell heterogeneity (het_seed; seed 0 = the production
cells), the L6 kick jitter, and the tone's fibre spike trains.

Each run is scored with the Stage 0 spike-based measure (volley_stats.py
--trains): cycles, duration and frequency of the train the kick evokes, TC /
RE participation per cycle, RE first-spike SD per volley.

    cd neuron
    ../.venv-neuron/bin/python thal_ring_test.py --g-gap 0.03,0.003,0.001,0.0005,0 \\
        --kick l6,tone --seeds 5 --jobs 4 --outdir ../res/2026-09-24/ring
"""

import argparse
import itertools
import json
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

MN5_THAL = {"tc": 346, "re": 91}
CORTEX = ["l4e", "l4i", "l23e", "l23i", "l5e", "l5i", "l6e", "l6i"]
G_L6_TC = G_L6_RE = 0.03        # production g_l6_tc / g_l6_re


def run_one(g_gap, seed, kick, kick_t=1000.0, tstop=4000.0, kick_sd=3.0,
            outdir=None, extra=None):
    """Build, kick once, run; returns (summary dict, npz path)."""
    os.chdir(HERE)   # compiled mechanisms live here
    from neuron import h
    import ctx_thalamus_mpi as M
    import auditory_input as AI
    import volley_stats as V

    sizes = dict(MN5_THAL, **{p: 0 for p in CORTEX})
    kw = dict(extra or {})
    net = M.ParallelCorticoThalamicNet(sizes=sizes, g_gap=g_gap, het_seed=seed, **kw)

    ncs, events = [], []
    if kick == "l6":
        for key, pop, gtot in [("l6_tc", "tc", G_L6_TC), ("l6_re", "re", G_L6_RE)]:
            lo, hi = net.ranges[pop]
            for gid in range(lo, hi):
                if (key, gid) not in net.syns:
                    continue
                rng = np.random.default_rng([seed, 5, gid])
                nc = h.NetCon(None, net.syns[(key, gid)])
                nc.weight[0] = gtot / 100.0
                ncs.append(nc)
                events.append((nc, kick_t + rng.normal(0.0, kick_sd, 100)))
        fih = h.FInitializeHandler(
            lambda: [nc.event(float(t)) for nc, ts in events for t in ts])
        aud = None
    elif kick == "tone":
        spec = {"seed": seed, "events": [{"t": kick_t, "type": "tone",
                                          "dur": 50.0, "level": 60.0}]}
        aud = AI.AuditoryInput(net, spec, tstop)
        fih = None
    else:
        raise ValueError(kick)

    wall = net.run(tstop=tstop)
    t, g = net.gather()
    t, g = np.asarray(t), np.asarray(g).astype(int)

    R = net.ranges
    trains = V.volley_trains(V.re_volleys(t, g, R, tstop, skip_ms=0.0))
    post = [tr for tr in trains if tr["t0"] >= kick_t]
    ev = post[0] if post and post[0]["t0"] < kick_t + 300.0 else None
    pre_spikes = int((t < kick_t).sum())
    s = {"g_gap": g_gap, "seed": seed, "kick": kick, "wall": wall,
         "params": dict(extra or {}),
         "pre_kick_spikes": pre_spikes,
         "tc_spikes": int(((g >= R["tc"][0]) & (g < R["tc"][1])).sum()),
         "re_spikes": int(((g >= R["re"][0]) & (g < R["re"][1])).sum()),
         "n_trains_after_kick": len(post),
         "evoked": ev is not None}
    if ev is not None:
        s.update({"cycles": ev["cycles"], "duration": ev["duration"],
                  "freq": ev["freq"], "wax_wane": ev["wax_wane"],
                  "re_part": [round(x, 3) for x in ev["re_part"]],
                  "tc_part": [round(x, 3) for x in ev["tc_part"]],
                  "re_sd_first": ev["re_sd"][0]})
    path = None
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        tag = "".join("_%s%g" % (k, v) for k, v in sorted((extra or {}).items()))
        path = os.path.join(outdir, "ring_%s_ggap%g%s_seed%d.npz" % (kick, g_gap, tag, seed))
        np.savez_compressed(path, times=t, gids=g, ranges=R, sizes=net.sizes,
                            tstop=tstop, kick_t=kick_t, kick=kick, g_gap=g_gap,
                            seed=seed, summary=json.dumps(s))
    net.teardown()
    return s, path


def _task(args):
    return run_one(*args)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--g-gap", default="0.03,0.003,0.001,0.0005,0")
    ap.add_argument("--kick", default="l6,tone")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--kick-t", type=float, default=1000.0)
    ap.add_argument("--tstop", type=float, default=4000.0)
    ap.add_argument("--kick-sd", type=float, default=3.0)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--outdir", default="")
    ap.add_argument("--param", action="append", default=[],
                    help="extra ParallelCorticoThalamicNet constructor grid, "
                         "repeatable: --param tau2_re_tc=8,20,40 --param g_tc_re=0.011,0.03")
    a = ap.parse_args()

    pnames = [p.split("=", 1)[0] for p in a.param]
    pvals = [[float(x) for x in p.split("=", 1)[1].split(",")] for p in a.param]
    combos = [dict(zip(pnames, c)) for c in itertools.product(*pvals)] if pnames else [{}]
    grid = list(itertools.product([k for k in a.kick.split(",")],
                                  [float(x) for x in a.g_gap.split(",")],
                                  combos, range(a.seeds)))
    tasks = [(gg, sd, k, a.kick_t, a.tstop, a.kick_sd, a.outdir or None, ex)
             for k, gg, ex, sd in grid]
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=a.jobs, mp_context=ctx,
                             max_tasks_per_child=1) as ex:
        results = [r for r, _ in ex.map(_task, tasks)]

    ptag = lambda s: " ".join("%s=%g" % kv for kv in sorted(s["params"].items()))
    print("%-5s %-7s %-4s | %-7s %-6s %-8s %-7s %-5s | %-9s %-22s | %s"
          % ("kick", "g_gap", "seed", "evoked", "cycles", "dur_ms", "Hz", "wax",
             "RE_SD_ms", "TC part per cycle", "trains after kick / pre-kick spikes"))
    print("-" * 118)
    for s in results:
        if s["evoked"]:
            tcp = ",".join("%.2f" % x for x in s["tc_part"][:6])
            print("%-5s %-7g %-4d | %-7s %-6d %-8.0f %-7s %-5s | %-9.2f %-22s | %d / %d"
                  % (s["kick"], s["g_gap"], s["seed"], "yes", s["cycles"], s["duration"],
                     "%.1f" % s["freq"] if s["cycles"] > 1 else "-",
                     "yes" if s["wax_wane"] else "no", s["re_sd_first"], tcp,
                     s["n_trains_after_kick"], s["pre_kick_spikes"]), ptag(s))
        else:
            print("%-5s %-7g %-4d | %-7s %-6s %-8s %-7s %-5s | %-9s %-22s | %d / %d"
                  % (s["kick"], s["g_gap"], s["seed"], "no", "-", "-", "-", "-", "-", "-",
                     s["n_trains_after_kick"], s["pre_kick_spikes"]), ptag(s))
    print("\nper (kick, g_gap), over seeds: evoked / median cycles / max cycles / "
          "runs with >= 6 cycles / median RE first-spike SD")
    for k, gg, ex in itertools.product(a.kick.split(","), [float(x) for x in a.g_gap.split(",")],
                                       combos):
        rs = [s for s in results if s["kick"] == k and s["g_gap"] == gg and s["params"] == ex]
        ev = [s for s in rs if s["evoked"]]
        cyc = [s["cycles"] for s in ev]
        print("  %-5s g_gap=%-7g %-28s evoked %d/%d, cycles median %s max %s, >=6: %d/%d, RE SD %s ms"
              % (k, gg, ptag(rs[0]) if rs else "", len(ev), len(rs), "%.0f" % np.median(cyc) if cyc else "-",
                 max(cyc) if cyc else "-", sum(c >= 6 for c in cyc), len(rs),
                 "%.2f" % np.median([s["re_sd_first"] for s in ev]) if ev else "-"))
    if a.outdir:
        with open(os.path.join(a.outdir, "ring_summary.json"), "w") as f:
            json.dump(results, f, indent=1, default=float)


if __name__ == "__main__":
    main()
