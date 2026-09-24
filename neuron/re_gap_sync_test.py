"""
Why is the RE population so synchronised? Isolates the RE<->RE gap junctions.

Builds the MN5-scale RE nucleus (91 tc_neuron.RECell, the same small-world
gap topology as ctx_thalamus_mpi._wire_gap_range: ring +-3 plus 2 random
shortcuts per cell, same seeds) WITHOUT TC, RE->RE inhibition or cortex, and
reports:

  1. a single RE cell's input conductance at rest, vs. one gap junction
  2. coupling coefficients (nearest neighbour / far side of the ring) for a
     range of g_gap values -- steady-state dV spread from a current step
  3. RE first-spike spread when every cell gets 100 inputs from a 693-cell
     "L6E" volley whose timing is deliberately jittered (SD sigma), with gap
     junctions on (g_gap) vs. off

Run from the repo root (needs the compiled mod files in neuron/):
    cd neuron && ../.venv-neuron/bin/python re_gap_sync_test.py
"""

import numpy as np
from neuron import h

import tc_neuron as T

h.load_file("stdrun.hoc")
h.celsius = 36.0
h.dt = 0.025

N_RE, N_L6, K, G_L6 = 91, 693, 100, 0.03     # MN5 scale 1.65, conv 100, g_l6_re
E_PAS_RE, HET = -82.0, 0.05


def _re_cells(het=True):
    cells = []
    for i in range(N_RE):
        c = T.RECell(gsk=1e-3)
        c.soma.e_pas = E_PAS_RE
        if het:  # same deterministic per-gid jitter scheme as ctx_thalamus_mpi._jitter
            c.soma.e_pas *= 1 + np.random.default_rng(10_000_000 + i * 97 + 1).uniform(-HET, HET)
            c.soma.gcabar_it2 *= 1 + np.random.default_rng(10_000_000 + i * 97 + 2).uniform(-HET, HET)
        cells.append(c)
    return cells


def _wire_gaps(cells, g, gap_deg=6, gap_short=2, seed_base=2000003):
    gaps, n, half = [], len(cells), gap_deg // 2
    if g <= 0:
        return gaps
    for i in range(n):
        nb = {(i + d) % n for d in range(1, half + 1)} | {(i - d) % n for d in range(1, half + 1)}
        r = np.random.default_rng(seed_base + i)
        nb |= {int(r.integers(0, n)) for _ in range(gap_short)}
        nb.discard(i)
        for j in nb:
            gp = h.Gap(cells[i].soma(0.5))
            gp.g = g
            h.setpointer(cells[j].soma(0.5)._ref_v, "vgap", gp)
            gaps.append(gp)
    return gaps


def _step_dv(cells, probe_idx, amp=-0.01, delay=2000.0, dur=1000.0):
    st = h.IClamp(cells[0].soma(0.5))
    st.delay, st.dur, st.amp = delay, dur, amp
    vs = [h.Vector().record(cells[k].soma(0.5)._ref_v) for k in probe_idx]
    tv = h.Vector().record(h._ref_t)
    h.finitialize(E_PAS_RE)
    h.continuerun(delay + dur)
    tv = np.array(tv)
    pre = (tv > delay - 100) & (tv < delay)
    post = (tv > delay + dur - 100) & (tv < delay + dur)
    return [np.array(v)[post].mean() - np.array(v)[pre].mean() for v in vs]


def input_conductance():
    c = T.RECell(gsk=1e-3)
    c.soma.e_pas = E_PAS_RE
    dv = _step_dv([c], [0])[0]
    return 0.01 / abs(dv)


def coupling(g):
    cells = _re_cells(het=False)
    _gaps = _wire_gaps(cells, g)
    d = _step_dv(cells, [0, 1, N_RE // 2])
    return d[1] / d[0], d[2] / d[0], len(_gaps) / N_RE


def volley_spread(sigma_ms, g, seed=1):
    cells = _re_cells()
    _gaps = _wire_gaps(cells, g)
    src_t = 1000.0 + np.random.default_rng(seed).normal(0, sigma_ms, N_L6)
    stims = []
    for tt in src_t:
        ns = h.NetStim()
        ns.number, ns.start, ns.noise = 1, max(1.0, tt), 0
        stims.append(ns)
    syns, ncs, recs = [], [], []
    for i, c in enumerate(cells):
        s = h.Exp2Syn(c.soma(0.5))
        s.e, s.tau1, s.tau2 = 0.0, 0.5, 2.0
        syns.append(s)
        for src in np.random.default_rng(i + 32).choice(N_L6, K, replace=False):
            nc = h.NetCon(stims[src], s)
            nc.weight[0], nc.delay = G_L6 / K, 1.0
            ncs.append(nc)
        v = h.Vector()
        det = h.NetCon(c.soma(0.5)._ref_v, None, sec=c.soma)
        det.threshold = -10
        det.record(v)
        recs.append((v, det))
    h.finitialize(E_PAS_RE)
    h.continuerun(1300)
    first = [np.array(v)[np.array(v) > 900].min() for v, _ in recs if np.any(np.array(v) > 900)]
    first = np.array(first)
    return len(first), (first.std() if len(first) > 1 else float("nan"))


def main():
    gin = input_conductance()
    print(f"RE input conductance at rest: {gin:.4f} uS; one gap junction (0.03 uS) = {0.03 / gin:.1f}x that")

    print("\ng_gap (uS) | junctions/cell | coupling: nearest neighbour | far side of ring")
    for g in (0.03, 0.01, 0.003, 0.001, 0.0005):
        a, b, per = coupling(g)
        print(f"  {g:7.4f}  |  {per:4.1f}          |  {a:.3f}                     |  {b:.3f}")

    print("\nL6 input spread (SD) | RE spread, gaps on (0.03 uS) | RE spread, gaps off")
    for sigma in (0.5, 3.0, 10.0, 25.0):
        n_on, s_on = volley_spread(sigma, 0.03)
        n_off, s_off = volley_spread(sigma, 0.0)
        print(f"  {sigma:5.1f} ms           |  {n_on:3d}/91 fire, {s_on:.3f} ms      "
              f"|  {n_off:3d}/91 fire, {s_off:.3f} ms")


if __name__ == "__main__":
    main()
