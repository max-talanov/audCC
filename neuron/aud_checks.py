"""
Stage A unit checks for the auditory input and the brain-state leak split
(PLAN-auditory-input.md, Stage A). No MPI; seconds on a laptop.

    ../.venv-neuron/bin/python aud_checks.py            # all checks
    ../.venv-neuron/bin/python aud_checks.py --only tc  # one check

Checks:
  fibres  -- fibre PSTH over many trials matches the rate model
  determ  -- fibre trains depend only on (seed, fibre id), not call order
  leak    -- nrem leak split reproduces the legacy voltage trace per cell type
  tc      -- driver calibration: coincident fibres needed to fire a TC cell in
             wake vs NREM, and burst vs tonic response to one tone
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from neuron import h                                    # noqa: E402

import auditory_input as AI                              # noqa: E402
import brain_state as BS                                 # noqa: E402
import tc_neuron as T                                    # noqa: E402
import cortex_neuron as C                                # noqa: E402

h.load_file("stdrun.hoc")


def check_fibres(n_trials=400, n_fibres=20):
    spec = AI.load_spec({"events": [{"t": 100.0 + 200.0 * i, "type": "tone",
                                     "dur": 50.0, "level": 60.0}
                                    for i in range(n_trials)]})
    tstop = 100.0 + 200.0 * n_trials
    s_edges = np.arange(-20.0, 101.0, 2.0)
    counts = np.zeros(len(s_edges) - 1)
    for fid in range(n_fibres):
        t = AI.fibre_train(fid, spec, tstop)
        rel = (t - 100.0) % 200.0
        rel[rel > 150] -= 200.0
        counts += np.histogram(rel, s_edges)[0]
    emp = counts / (n_trials * n_fibres) / (np.diff(s_edges) / 1000.0)
    mid = 0.5 * (s_edges[:-1] + s_edges[1:])
    fine = np.arange(-20.0, 100.0, 0.1)
    model_fine = AI.drive_amplitude(60.0, spec) * AI.rate_envelope(fine, 50.0, spec)
    model = np.array([model_fine[(fine >= a) & (fine < b)].mean()
                      for a, b in zip(s_edges[:-1], s_edges[1:])])
    r = np.corrcoef(emp, model)[0, 1]
    ratio = emp.sum() / model.sum()
    print("[fibres] peak model %.0f Hz, empirical %.0f Hz; corr %.3f; "
          "total empirical/model %.3f (refractoriness lowers it slightly)"
          % (model.max(), emp.max(), r, ratio))
    ok = r > 0.98 and 0.85 < ratio < 1.02
    print("[fibres] %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_determinism():
    spec = AI.load_spec({"protocol": {"type": "random_tones", "start": 1000,
                                      "stop": 20000, "isi_min": 2000,
                                      "isi_max": 4000, "dur": 50, "level": 60},
                         "r_spont": 3.0})
    a = [AI.fibre_train(f, spec, 20000.0) for f in range(10)]
    b = [AI.fibre_train(f, spec, 20000.0) for f in reversed(range(10))][::-1]
    same_trains = all(np.array_equal(x, y) for x, y in zip(a, b))
    same_map = np.array_equal(AI.fibres_of_tc(5, spec, 40), AI.fibres_of_tc(5, spec, 40))
    spec2 = AI.load_spec({"protocol": {"type": "random_tones", "start": 1000,
                                       "stop": 20000, "isi_min": 2000,
                                       "isi_max": 4000, "dur": 50, "level": 60},
                          "r_spont": 3.0})
    same_sched = [e["t"] for e in spec["events"]] == [e["t"] for e in spec2["events"]]
    ok = same_trains and same_map and same_sched
    print("[determ] trains order-independent: %s; TC->fibre map stable: %s; "
          "schedule reproducible: %s (%d tones) -> %s"
          % (same_trains, same_map, same_sched, len(spec["events"]),
             "PASS" if ok else "FAIL"))
    return ok


def _make(kind, e_pas):
    if kind == "tc":
        c = T.TCCell(gsk=0.0, gh=0.0)
        c.soma.e_pas = e_pas
    elif kind == "re":
        c = T.RECell(gsk=1e-3)
        c.soma.e_pas = e_pas
    elif kind == "ib":
        c = C.PYCellIB(e_pas=e_pas, gnap=2e-4)
    elif kind == "py":
        c = C.PYCell(e_pas=e_pas, gsk=8e-4)
    else:
        c = C.FSCell(e_pas=e_pas)
    return c


def _trace(kind, e_pas, state, tstop=1500.0, amp=None):
    c = _make(kind, e_pas)
    pop = {"tc": "tc", "re": "re", "ib": "l5e", "py": "l4e", "fs": "l4i"}[kind]
    if state:
        BS.apply_state(c, pop, state)
    # a hyperpolarising then depolarising step exercises the leak, I_T
    # rebound, and spiking, so any mismatch in the split shows up.
    stims = []
    for t0, a in [(300.0, -0.1 if amp is None else amp), (900.0, 0.15)]:
        ic = h.IClamp(c.soma(0.5))
        ic.delay, ic.dur, ic.amp = t0, 200.0, a
        stims.append(ic)
    v = h.Vector().record(c.soma(0.5)._ref_v)
    h.dt, h.celsius = 0.025, 36.0
    h.finitialize(-70)
    h.continuerun(tstop)
    return np.array(v)


def check_leak():
    ok = True
    for kind, e_pas in [("tc", -80.0), ("re", -82.0), ("ib", -70.0),
                        ("py", -70.0), ("fs", -68.0)]:
        v0 = _trace(kind, e_pas, None)
        v1 = _trace(kind, e_pas, "nrem")
        d = np.abs(v0 - v1).max()
        nsp0 = int(((v0[1:] > -10) & (v0[:-1] <= -10)).sum())
        nsp1 = int(((v1[1:] > -10) & (v1[:-1] <= -10)).sum())
        good = d < 1e-6 and nsp0 == nsp1
        ok &= good
        print("[leak] %-3s max |dV| legacy vs nrem = %.2e mV, spikes %d vs %d  %s"
              % (kind, d, nsp0, nsp1, "PASS" if good else "FAIL"))
    return ok


def _rest(kind, e_pas, state, t=2000.0):
    h.dt, h.celsius = 0.025, 36.0
    c = _make(kind, e_pas)
    pop = {"tc": "tc", "re": "re", "ib": "l5e", "py": "l4e", "fs": "l4i"}[kind]
    BS.apply_state(c, pop, state)
    v = h.Vector().record(c.soma(0.5)._ref_v)
    h.finitialize(-70)
    h.continuerun(t)
    return float(np.array(v)[-1])


def _tc_cell(state, hold=None):
    """Network TC cell (e_pas -80, leak split for `state`), optionally held
    at `hold` mV by a DC current found by bisection."""
    h.dt, h.celsius = 0.025, 36.0
    c = _make("tc", -80.0)
    BS.apply_state(c, "tc", state)
    c.dc = None
    if hold is not None:
        c.dc = h.IClamp(c.soma(0.5))
        c.dc.delay, c.dc.dur = 0.0, 1e9
        lo, hi = -1.0, 1.0
        for _ in range(16):
            c.dc.amp = 0.5 * (lo + hi)
            h.dt, h.celsius = 0.025, 36.0
            h.finitialize(-70)
            h.continuerun(800.0)
            if c.soma(0.5).v > hold:
                hi = c.dc.amp
            else:
                lo = c.dc.amp
    return c


def _tc_response(state, w, n_coinc=None, tone=False, hold=None, na=True,
                 tstop=1300.0):
    """Drive a TC cell at 1000 ms with n_coinc coincident fibre spikes, or
    with one 60 dB tone through k fibres. Returns (spike times, soma V)."""
    c = _tc_cell(state, hold)
    if not na:
        c.soma.gnabar_hh2 = 0.0
    syn = h.Exp2Syn(c.soma(0.5))
    syn.e, syn.tau1, syn.tau2 = 0.0, AI.DEFAULTS["tau1"], AI.DEFAULTS["tau2"]
    if tone:
        spec = AI.load_spec({"events": [{"t": 1000.0, "type": "tone",
                                         "dur": 50.0, "level": 60.0}], "w": w})
        trains = [AI.fibre_train(f, spec, tstop) for f in AI.fibres_of_tc(0, spec, 50)]
    else:
        trains = [np.array([1000.0])] * n_coinc
    ncs = []
    for tr in trains:
        nc = h.NetCon(None, syn)
        nc.weight[0] = w
        ncs.append((nc, tr))
    fih = h.FInitializeHandler(lambda: [nc.event(float(t)) for nc, tr in ncs for t in tr])
    spk = h.Vector()
    det = h.NetCon(c.soma(0.5)._ref_v, None, sec=c.soma)
    det.threshold = -10
    det.record(spk)
    v = h.Vector().record(c.soma(0.5)._ref_v)
    tv = h.Vector().record(h._ref_t)
    h.dt, h.celsius = 0.025, 36.0
    h.finitialize(-70)
    h.continuerun(tstop)
    del fih
    return np.array(spk), np.array(tv), np.array(v)


def _is_burst(spk, t0=1000.0):
    """Lu et al. 1992: first evoked spike preceded by >=100 ms silence and
    followed by an ISI <= 4 ms."""
    ev = spk[spk >= t0]
    if len(ev) < 2:
        return False
    prev = spk[spk < ev[0]]
    silent = len(prev) == 0 or ev[0] - prev[-1] >= 100.0
    return silent and (ev[1] - ev[0]) <= 4.0


def check_tc(ws=(0.001, 0.002, 0.003, 0.004, 0.006), holds=(-68.0, -80.0)):
    # "Depolarised" is -68 mV, not the -60 mV of PLAN Stage A: this TCCell
    # fires tonically on its own from about -67 mV (DC sweep, e_pas -80),
    # so it cannot be held higher without spontaneous spiking.
    for st in ("nrem", "wake"):
        spk, tv, v = _tc_response(st, 0.0, n_coinc=0, tstop=3000.0)
        n = int((spk >= 1000.0).sum())
        print("[tc] %-4s state, no input: %s"
              % (st, "%.1f Hz spontaneous firing (not at rest)" % (n / 2.0) if n
                 else "silent, rest %.1f mV" % v[-1]))
    print("[tc] single-fibre EPSP (Na blocked) and spikes evoked by n coincident")
    print("     fibre spikes in 0-50 ms, cell DC-held at %s mV" % "/".join("%g" % x for x in holds))
    hdr = " | ".join("EPSP  n=1/2/3/4 @%g" % x for x in holds)
    print("     %-7s | %s" % ("w (uS)", hdr))
    for w in ws:
        cols = []
        for hold in holds:
            _, tv, v = _tc_response("nrem", w, 1, hold=hold, na=False)
            m = (tv >= 1000) & (tv < 1050)
            epsp = v[m].max() - v[(tv >= 990) & (tv < 1000)].mean()
            n = "/".join(str(int(((s >= 1000) & (s < 1050)).sum()))
                         for s in (_tc_response("nrem", w, k, hold=hold)[0]
                                   for k in (1, 2, 3, 4)))
            cols.append("%4.1f mV  %-9s" % (epsp, n))
        print("     %-7.4f | %s" % (w, " | ".join(cols)))
    w = AI.DEFAULTS["w"]
    print("[tc] one 60 dB tone via k=%d fibres, w=%.4f:" % (AI.DEFAULTS["k"], w))
    out = {}
    for label, st, hold in [("held -68", "nrem", -68.0), ("held -80", "nrem", -80.0),
                            ("wake state", "wake", None), ("nrem state", "nrem", None)]:
        s = _tc_response(st, w, tone=True, hold=hold)[0]
        ev = s[(s >= 1000) & (s < 1100)]
        isi = np.diff(ev)
        out[label] = (len(ev), _is_burst(s))
        if st == "wake" and hold is None:
            label += "*"
        print("     %-10s: %2d spikes in 0-100 ms, first at %5s ms, min ISI %5s ms, "
              "burst (Lu 1992): %s"
              % (label, len(ev), "%.1f" % (ev[0] - 1000) if len(ev) else "-",
                 "%.1f" % isi.min() if len(isi) else "-", _is_burst(s)))
    print("     * the placeholder wake preset fires spontaneously (see above), "
          "so this row is not a pure tone response")
    ok = (not out["held -68"][1]) and out["held -80"][1]
    print("[tc] tonic at -68 mV and burst at -80 mV: %s" % ("PASS" if ok else "FAIL"))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="fibres, determ, leak or tc")
    a = ap.parse_args()
    checks = {"fibres": check_fibres, "determ": check_determinism,
              "leak": check_leak, "tc": check_tc}
    for name, fn in checks.items():
        if not a.only or a.only == name:
            fn()


if __name__ == "__main__":
    main()
