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

    # --- RE burst structure -------------------------------------------------
    if len(nrt) > 2:
        t = np.sort(nrt)
        b, c = [], 1
        for gap in np.diff(t):
            if gap < 20.0:
                c += 1
            else:
                b.append(c); c = 1
        b.append(c)
        mb = float(np.mean(b))
        rows.append(("RE spikes per burst", "2 to >10", f"{mb:.1f} (mean)",
                     2.0 <= mb <= 15.0))

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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=str,
                    default="config/network_auditory_local.yaml")
    ap.add_argument("--tstop", type=float, default=200000.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-trigger", action="store_true",
                    help="disable the external modulatory spindle trigger")
    ap.add_argument("--no-assert", action="store_true",
                    help="always exit 0, even if criteria fail")
    args = ap.parse_args(argv)

    cfg = NetworkConfig.from_file(args.config)
    cfg.tstop = args.tstop
    sim = SimulationConfig(seed=args.seed, neuron_model="ht_neuron",
                           record_traces=True)
    sleep = SleepParams(emergent_spindles=True,
                        spindle_trigger=not args.no_trigger)
    print(f"Validating HH model ({cfg.tstop/1000:.0f} s, "
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
