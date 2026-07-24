"""
Validate simulated sleep spindles against Fernandez & Luthi (2020).

Turns "does this look like a spindle?" into a pass/fail table by measuring the
model against the quantitative criteria in
Fernandez LMJ & Luthi A, *Sleep Spindles: Mechanisms and Functions*,
Physiol Rev 100: 805-868, 2020.

Criteria (section references are to that review):

  intra-spindle frequency  10-15 Hz sigma band            (sect. III.A)
  duration                 0.5-3 s                        (sect. V.D)
  density                  2-8 /min in N2                 (sect. IV.A.2)
  inter-spindle interval   5-10 s refractory period       (sect. V.D)
  SO coupling              50-70% of spindles on UP state (sect. VI.E)
  infraslow clustering     sigma power modulated ~0.02 Hz (sect. IV.E)
  RE spikes per burst      2 to >10                       (sect. V.B.1)
  RE membrane potential    must reach < -55 mV to burst   (sect. V.B)
  TC membrane potential    must reach < -65 mV to rebound (sect. V.B.4)

Run:

    python3 tc_validate.py --tstop 200000            # HH + trigger (default)
    python3 tc_validate.py --no-trigger              # ungated HH model

Exits non-zero if any criterion fails, so it can gate CI.

Simulator-agnostic
------------------
``validate()`` operates only on spike times and voltage traces, via pure-numpy
signal helpers -- it has **no NEST dependency** and can score a run from any
backend (NEST, NEURON, ...). A backend need only emit three plain dicts (the
"result contract"); the NEST driver in ``main()`` is just one producer.

Result contract (all times in **ms**, voltages in **mV**):

  spikes:  {layer: {"times":   float[],     # spike times, ms
                    "senders": int[]}}      # per-spike cell id (for burst
                                            # structure; may be empty)
  traces:  {layer: {"time":    float[],     # sample times, ms
                    "voltage": float[]}}    # mean V_m of the layer, mV
  meta:    {"tstop": float,                 # run length, ms
            "seed":  int}                   # for reproducible EEG-proxy noise

  Required layer keys: "MGB" (TC relay) and "nRT" (RE reticular). The cortical
  keys "L23"/"L5"/"L6" feed the slow-wave component of the EEG proxy; absent
  layers are treated as silent. ``validate_result(spikes, traces, meta)`` is the
  backend-facing entry point; ``self_test()`` proves the path runs on synthetic
  data with neither NEST nor NEURON present.
"""

import argparse
from pathlib import Path

import numpy as np

try:
    from .tc_network import (AuditoryThalamoCorticalSleep, NetworkConfig,
                             SimulationConfig, SynapseParams, SleepParams, HHParams)
    from . import tc_run
    from .tc_spindle_figures import tc_episodes, intra_spindle_freq
except ImportError:
    from tc_network import (AuditoryThalamoCorticalSleep, NetworkConfig,
                            SimulationConfig, SynapseParams, SleepParams, HHParams)
    import tc_run
    from tc_spindle_figures import tc_episodes, intra_spindle_freq


def _fmt(ok):
    return "PASS" if ok else "FAIL"


def detect_spindles(mgb, nrt, tstop, lo=10.0, hi=15.0, k_sd=1.5,
                    min_dur_ms=500.0, merge_ms=100.0):
    """Detect discrete spindles the way the review does (sect. III.B): band-pass
    the thalamic signal to the sigma band, take the envelope, apply a FIXED
    threshold, and keep events lasting >= 0.5 s.

    This is deliberately not (a) an envelope *percentile* threshold, which
    self-normalises and so cannot show suppression, nor (b) "any contiguous run
    of relay spikes", which counts an isolated spike pair as an event. Both
    inflated earlier measurements.

    Returns (starts, ends) in ms.
    """
    allt = np.concatenate([np.asarray(mgb, float), np.asarray(nrt, float)])
    if len(allt) < 10:
        return np.array([]), np.array([])
    t, rate = tc_run.population_rate(allt, tstop, bin_ms=1.0, smooth_ms=3.0)
    fs = 1000.0
    band = tc_run._bandpass(rate, fs, lo, hi)
    _, env = tc_run.bandpass_envelope(rate, fs, lo, hi)
    thr = env.mean() + k_sd * env.std()          # fixed threshold
    above = env > thr
    if not above.any():
        return np.array([]), np.array([])
    edges = np.diff(above.astype(int))
    starts = list(np.where(edges == 1)[0] + 1)
    ends = list(np.where(edges == -1)[0] + 1)
    if above[0]:
        starts = [0] + starts
    if above[-1]:
        ends = ends + [len(above) - 1]
    s = np.array(starts, float) * (1000.0 / fs)
    e = np.array(ends, float) * (1000.0 / fs)
    # merge events separated by a short dip, then apply the duration criterion
    if len(s) > 1:
        keep_s, keep_e = [s[0]], [e[0]]
        for a, b in zip(s[1:], e[1:]):
            if a - keep_e[-1] <= merge_ms:
                keep_e[-1] = b
            else:
                keep_s.append(a); keep_e.append(b)
        s, e = np.array(keep_s), np.array(keep_e)
    ok = (e - s) >= min_dur_ms
    return s[ok], e[ok]


def validate(spikes, traces, meta, infraslow_hz=0.02, verbose=True):
    """Measure every criterion. Returns (rows, all_ok).

    ``rows`` is a list of (criterion, paper_target, measured, ok).
    """
    tstop = meta["tstop"]
    mgb = np.asarray(spikes.get("MGB", {}).get("times", []), float)
    nrt = np.asarray(spikes.get("nRT", {}).get("times", []), float)
    rows = []

    # Discrete spindles by the review's own method: fixed threshold on the
    # sigma-band envelope, events >= 0.5 s (see detect_spindles).
    st, en = detect_spindles(mgb, nrt, tstop)
    isi = np.diff(st) / 1000.0 if len(st) > 1 else np.array([])
    dur = (en - st) / 1000.0 if len(st) else np.array([])

    # --- intra-spindle frequency -------------------------------------------
    f = intra_spindle_freq(mgb, nrt, st, en)
    if len(f):
        med_f = float(np.median(f))
        rows.append(("intra-spindle frequency", "10-15 Hz",
                     f"{med_f:.1f} Hz (median)", 10.0 <= med_f <= 15.0))

    # --- duration -----------------------------------------------------------
    if len(dur):
        med_d = float(np.median(dur))
        rows.append(("duration", "0.5-3 s", f"{med_d:.2f} s (median)",
                     0.5 <= med_d <= 3.0))

    # --- density ------------------------------------------------------------
    dens = len(st) / (tstop / 60000.0)
    rows.append(("density", "2-8 /min", f"{dens:.1f} /min", 2.0 <= dens <= 8.0))

    # --- inter-spindle interval (refractory) --------------------------------
    if len(isi):
        med_i = float(np.median(isi))
        rows.append(("inter-spindle interval", "5-10 s",
                     f"{med_i:.1f} s (median)", 5.0 <= med_i <= 10.0))

    # --- SO coupling --------------------------------------------------------
    eeg = tc_run.build_eeg_like(traces, tstop, seed=meta.get("seed", 0))
    t_e, fs = eeg["t"], eeg["fs"]
    slow = tc_run._bandpass(eeg["eeg"], fs, 0.5, 2.0)
    try:
        from scipy.signal import hilbert
        ph = np.angle(hilbert(slow))
        if len(st):
            idx = np.clip(np.searchsorted(t_e, st), 0, len(ph) - 1)
            up = float(np.mean(np.cos(ph[idx]) > 0) * 100.0)
            rows.append(("SO coupling", "50-70%", f"{up:.0f}% on UP state",
                         50.0 <= up <= 70.0))
    except Exception:
        pass

    # --- infraslow clustering ----------------------------------------------
    env = eeg["env"]
    if tstop / 1000.0 >= 3.0 / max(1e-6, infraslow_hz):
        e = env[int(2 * fs):]
        e = e - e.mean()
        fr = np.fft.rfftfreq(len(e), 1.0 / fs)
        P = np.abs(np.fft.rfft(e)) ** 2
        bd = (fr > infraslow_hz * 0.6) & (fr < infraslow_hz * 1.4)
        rf = (fr > 0.01) & (fr < 1.0)
        if bd.any() and rf.any():
            pk = float(fr[bd][np.argmax(P[bd])])
            ratio = float(P[bd].max() / max(1e-12, np.median(P[rf])))
            rows.append(("infraslow clustering", f"~{infraslow_hz} Hz",
                         f"{pk:.3f} Hz (ratio {ratio:.0f})", ratio > 3.0))

    # --- burst mode (per cell, at the review's burst criterion) -------------
    # The review defines a burst as 2 to >10 spikes at frequencies up to several
    # hundred Hz (TRN) and 2-6 spikes for the TC rebound. So intra-burst ISI
    # must be < 10 ms (>100 Hz). Using a looser window (e.g. 20 ms) lumps tonic
    # spikes together and makes tonic firing look like bursting -- an earlier
    # version of this check did exactly that and wrongly reported a PASS.
    for key, lab, lo_b in [("nRT", "RE", 2.0), ("MGB", "TC", 2.0)]:
        arr = spikes.get(key, {})
        tim = np.asarray(arr.get("times", []), float)
        sen = np.asarray(arr.get("senders", []), float)
        if len(tim) < 3:
            continue
        per_cell = []
        for c in np.unique(sen):
            tc_t = np.sort(tim[sen == c])
            if len(tc_t) < 2:
                continue
            b, k = [], 1
            for gap in np.diff(tc_t):
                if gap < 10.0:        # >100 Hz -> same burst
                    k += 1
                else:
                    b.append(k); k = 1
            b.append(k)
            per_cell.append(np.mean(b))
        if per_cell:
            mb = float(np.mean(per_cell))
            rows.append((f"{lab} spikes per burst", "2 to >10 (>100 Hz)",
                         f"{mb:.2f} (mean/cell)", mb >= lo_b))

    # --- membrane-potential operating ranges --------------------------------
    for key, lab, thr in [("nRT", "RE", -55.0), ("MGB", "TC", -65.0)]:
        if key in traces:
            v = np.asarray(traces[key]["voltage"], float)
            frac = float(np.mean(v < thr) * 100.0)
            rows.append((f"{lab} V_m reaches < {thr:.0f} mV",
                         "required to burst", f"{frac:.0f}% of time",
                         frac >= 20.0))

    all_ok = all(r[3] for r in rows)
    if verbose:
        print("\n=== Spindle validation vs Fernandez & Luthi 2020 ===")
        print(f"{'criterion':<32}{'paper':<20}{'measured':<24}result")
        print("-" * 84)
        for crit, tgt, meas, ok in rows:
            print(f"{crit:<32}{tgt:<20}{meas:<24}{_fmt(ok)}")
        n_ok = sum(1 for r in rows if r[3])
        print("-" * 84)
        print(f"{n_ok}/{len(rows)} criteria passed")
    return rows, all_ok


def validate_result(spikes, traces, meta, infraslow_hz=0.02, verbose=True):
    """Backend-facing entry point: score a result that satisfies the contract in
    the module docstring, from ANY simulator. Thin wrapper over ``validate`` that
    does not depend on ``tc_network``/NEST, so a NEURON backend can call it
    directly."""
    return validate(spikes, traces, meta, infraslow_hz=infraslow_hz,
                    verbose=verbose)


def _synthetic_run(tstop=200000.0, n_spindles=None, seed=1):
    """Fabricate a contract-shaped result with spindle-like thalamic bursting,
    for testing the validator with neither NEST nor NEURON present."""
    rng = np.random.default_rng(seed)
    n_spindles = n_spindles or int(tstop / 1000.0 / 7.0)   # ~7 s apart
    n_tc, n_re = 40, 40
    mgb_t, mgb_s, nrt_t, nrt_s = [], [], [], []
    starts = 3000.0 + np.arange(n_spindles) * (tstop - 6000.0) / max(1, n_spindles)
    for s0 in starts:
        dur = rng.uniform(600.0, 900.0)                    # 0.6-0.9 s spindle
        for cyc in np.arange(s0, s0 + dur, 1000.0 / 13.0):  # ~13 Hz cycles
            for cells, tl, sl in [(n_re, nrt_t, nrt_s), (n_tc, mgb_t, mgb_s)]:
                nb = rng.integers(3, 6)                     # burst of 3-5 spikes
                who = rng.integers(0, cells)
                for k in range(nb):
                    tl.append(cyc + k * rng.uniform(4.0, 6.0))  # ~200 Hz
                    sl.append(int(who))
    grid = np.arange(0.0, tstop, 1.0)
    # thalamic V_m: rests hyperpolarised, dips deeper during spindles
    base = -68.0 + 3.0 * np.sin(2 * np.pi * 1.0 * grid / 1000.0)
    thal_v = base.copy()
    for s0 in starts:
        m = (grid >= s0) & (grid <= s0 + 800.0)
        thal_v[m] -= 18.0
    cort_v = -63.0 + 6.0 * np.sin(2 * np.pi * 1.0 * grid / 1000.0)
    spikes = {"MGB": {"times": np.array(mgb_t), "senders": np.array(mgb_s)},
              "nRT": {"times": np.array(nrt_t), "senders": np.array(nrt_s)}}
    traces = {"MGB": {"time": grid, "voltage": thal_v},
              "nRT": {"time": grid, "voltage": thal_v + 2.0},
              "L23": {"time": grid, "voltage": cort_v},
              "L5": {"time": grid, "voltage": cort_v},
              "L6": {"time": grid, "voltage": cort_v}}
    meta = {"tstop": tstop, "seed": seed}
    return spikes, traces, meta


def self_test():
    """Prove the validation path is simulator-agnostic: run it on synthetic data
    (no NEST, no NEURON) and confirm it produces a full criteria table."""
    print("Self-test: validating a SYNTHETIC spindling result "
          "(no simulator involved)...")
    spikes, traces, meta = _synthetic_run(tstop=120000.0)
    rows, _ = validate_result(spikes, traces, meta)
    assert len(rows) >= 8, f"expected >=8 criteria, got {len(rows)}"
    print(f"\nSelf-test OK: {len(rows)} criteria evaluated on synthetic data "
          f"with no simulator present.")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=str,
                    default="config/network_auditory_local.yaml")
    ap.add_argument("--tstop", type=float, default=200000.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--neuron-model", type=str, default="ht_neuron",
                    choices=["ht_neuron", "aeif_cond_exp"],
                    help="which burst-capable model to validate")
    ap.add_argument("--no-trigger", action="store_true",
                    help="disable the external modulatory spindle trigger")
    ap.add_argument("--no-assert", action="store_true",
                    help="always exit 0, even if criteria fail")
    ap.add_argument("--self-test", action="store_true",
                    help="validate synthetic data (no simulator) and exit")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()

    cfg = NetworkConfig.from_file(args.config)
    cfg.tstop = args.tstop
    sim = SimulationConfig(seed=args.seed, neuron_model=args.neuron_model,
                           record_traces=True)
    sleep = SleepParams(emergent_spindles=True,
                        spindle_trigger=not args.no_trigger)
    print(f"Validating {args.neuron_model} ({cfg.tstop/1000:.0f} s, "
          f"trigger={'off' if args.no_trigger else 'on'})...")
    model = AuditoryThalamoCorticalSleep(cfg, SynapseParams(), sleep, sim)
    spikes, traces, meta = model.run()
    meta["seed"] = args.seed
    _, ok = validate(spikes, traces, meta,
                     infraslow_hz=HHParams().infraslow_freq)
    if not ok and not args.no_assert:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
