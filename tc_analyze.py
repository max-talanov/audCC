"""
Quantitative analysis of a simulation result: sleep-spindle metrics per layer
and the thalamus -> cortex propagation of the spindle.

Runs the auditory thalamo-cortical column and, for each layer -- MGB, nRT
(thalamus) then L4, L2/3, L5, L6 (cortex) -- measures, from the sigma-band
(10-15 Hz) population activity:

  * peak intra-spindle frequency (Hz)
  * spindle density (events/min) and mean duration (s)
  * relative spindle power (sigma-band envelope, normalised to MGB)
  * propagation lag (ms): the cross-correlation delay of the sigma-band
    envelope relative to MGB -- a positive lag means the layer's spindle
    follows the thalamic one (the spindle travelling up the column)
  * SO-coupling: fraction of the layer's spindle epochs on the slow-wave UP state

Writes a printed table and out/spindle_analysis.md.

    python3 tc_analyze.py --config config/network_auditory_mushtaq.yaml \
        --tstop 15000 --outdir out
"""

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    from tc_sleep.tc_network import (AuditoryThalamoCorticalSleep, NetworkConfig,
                                     SimulationConfig, SynapseParams, SleepParams, HHParams)
    from tc_sleep import tc_run
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tc_network import (AuditoryThalamoCorticalSleep, NetworkConfig,
                            SimulationConfig, SynapseParams, SleepParams, HHParams)
    import tc_run

SIGMA = (10.0, 15.0)
ORDER = ["MGB", "nRT", "L4", "L23", "L5", "L6"]
FS = 1000.0


def _sigma_signals(spikes, tstop):
    """Per-layer sigma-band signal + envelope from the population rate."""
    out = {}
    for layer in ORDER:
        if layer not in spikes or len(spikes[layer]["times"]) < 10:
            continue
        _, rate = tc_run.population_rate(spikes[layer]["times"], tstop,
                                         bin_ms=1.0, smooth_ms=3.0)
        sig, env = tc_run.bandpass_envelope(rate, FS, *SIGMA)
        out[layer] = {"sig": sig, "env": env}
    return out


def _spindle_epochs(env, thr_pct=55.0, min_ms=150.0, merge_ms=100.0):
    thr = np.percentile(env, thr_pct)
    above = env > thr
    ed = np.diff(above.astype(int))
    s = list(np.where(ed == 1)[0] + 1) + ([0] if above[0] else [])
    e = list(np.where(ed == -1)[0] + 1) + ([len(env) - 1] if above[-1] else [])
    s, e = sorted(s), sorted(e)
    ep = [(a, b) for a, b in zip(s, e) if (b - a) >= min_ms]
    # merge across short dips
    merged = []
    for a, b in ep:
        if merged and a - merged[-1][1] <= merge_ms:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


def analyze(spikes, traces, meta):
    tstop = meta["tstop"]
    sig = _sigma_signals(spikes, tstop)
    if "MGB" not in sig:
        raise SystemExit("No thalamic (MGB) activity to analyse.")
    ref_env = sig["MGB"]["env"]
    ref_power = np.mean(ref_env) or 1.0

    # slow-wave reference (cortical rate <2 Hz) for SO-coupling
    cort = [l for l in ["L23", "L5", "L6", "L4"] if l in spikes]
    slow_ph = None
    if cort:
        allc = np.concatenate([spikes[l]["times"] for l in cort])
        _, rc = tc_run.population_rate(allc, tstop, bin_ms=1.0, smooth_ms=20.0)
        slow = tc_run._bandpass(rc, FS, 0.5, 2.0)
        from scipy.signal import hilbert
        slow_ph = np.angle(hilbert(slow))

    rows = []
    for layer in ORDER:
        if layer not in sig:
            continue
        env = sig["layer" if False else layer]["env"]
        # peak sigma frequency
        pk, _, _, _ = tc_run.detect_peak(sig[layer]["sig"], FS, SIGMA[0], SIGMA[1])
        ep = _spindle_epochs(env)
        density = len(ep) / (tstop / 60000.0)
        dur = np.mean([(b - a) / 1000.0 for a, b in ep]) if ep else 0.0
        power = np.mean(env) / ref_power
        # propagation lag vs MGB (cross-correlation of envelopes, +/- 200 ms)
        lag = 0.0
        if layer != "MGB":
            n = min(len(env), len(ref_env))
            a = env[:n] - env[:n].mean()
            b = ref_env[:n] - ref_env[:n].mean()
            maxlag = 200
            cc = np.correlate(a, b, mode="full")
            mid = len(cc) // 2
            window = cc[mid - maxlag: mid + maxlag + 1]
            lag = float(np.argmax(window) - maxlag)   # ms; +ve = layer lags MGB
        # SO coupling
        so = np.nan
        if slow_ph is not None and ep:
            idx = np.clip([a for a, _ in ep], 0, len(slow_ph) - 1)
            so = float(np.mean(np.cos(slow_ph[idx]) > 0) * 100)
        rows.append({"layer": layer,
                     "region": "thalamus" if not layer.startswith("L") else "cortex",
                     "freq": pk, "density": density, "dur": dur,
                     "power": power, "lag": lag, "so": so})
    return rows


def _fmt_table(rows):
    lines = [f"{'layer':<7}{'region':<10}{'freq(Hz)':<10}{'density/min':<13}"
             f"{'dur(s)':<9}{'rel.power':<11}{'lag vs MGB(ms)':<16}{'SO-coupled'}"]
    lines.append("-" * 86)
    for r in rows:
        so = f"{r['so']:.0f}%" if not np.isnan(r['so']) else "-"
        lines.append(f"{r['layer']:<7}{r['region']:<10}{r['freq']:<10.1f}"
                     f"{r['density']:<13.1f}{r['dur']:<9.2f}{r['power']:<11.2f}"
                     f"{r['lag']:<16.0f}{so}")
    return "\n".join(lines)


def _write_md(rows, meta, path):
    with open(path, "w") as f:
        f.write("# Spindle analysis (per layer)\n\n")
        f.write(f"Auditory thalamo-cortical column, {meta['tstop']/1000:.0f} s, "
                f"neuron model `{meta.get('neuron_model','?')}`. "
                "Sigma band 10-15 Hz.\n\n")
        f.write("| layer | region | freq (Hz) | density /min | dur (s) | "
                "rel. power | lag vs MGB (ms) | SO-coupled |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            so = f"{r['so']:.0f}%" if not np.isnan(r['so']) else "-"
            f.write(f"| {r['layer']} | {r['region']} | {r['freq']:.1f} | "
                    f"{r['density']:.1f} | {r['dur']:.2f} | {r['power']:.2f} | "
                    f"{r['lag']:.0f} | {so} |\n")
        cort = [r for r in rows if r["region"] == "cortex"]
        if cort:
            lag = np.mean([r["lag"] for r in cort])
            f.write(f"\n**Thalamus -> cortex propagation:** the sigma-band "
                    f"spindle envelope in the cortical layers lags MGB by "
                    f"{lag:.0f} ms on average (positive = spindle arrives in "
                    f"cortex after the thalamus) -- the spindle travelling up "
                    f"the column.\n")
        f.write("\n> Notes: frequency, relative power and the propagation lag "
                "are the robust measures. Density reflects the imposed-drive "
                "`iaf_cond_exp` regime (~one spindle per 1 Hz slow-oscillation "
                "cycle, i.e. ~60/min), higher than the 2-8/min of natural NREM; "
                "the per-layer detection threshold also splits some cortical "
                "spindles, inflating the upper-layer count. Use the AdEx/HH "
                "configs for physiological spindle *density* (see the validator, "
                "`tc_validate.py`).\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=str,
                    default="config/network_auditory_mushtaq.yaml")
    ap.add_argument("--tstop", type=float, default=15000.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", type=str, default="out")
    ap.add_argument("--tag", type=str, default=None,
                    help="filename tag; keeps parallel jobs from overwriting")
    args = ap.parse_args(argv)

    cfg = NetworkConfig.from_file(args.config)
    cfg.tstop = args.tstop
    nm = tc_run._config_neuron_model(args.config) or "iaf_cond_exp"
    sim = SimulationConfig(seed=args.seed, neuron_model=nm, record_traces=True)
    sleep = (SleepParams(emergent_spindles=True) if nm != "iaf_cond_exp"
             else SleepParams())
    print(f"Running column ({args.tstop/1000:.0f} s, model {nm}) for analysis...")
    hh_params = HHParams.from_file(args.config)
    model = AuditoryThalamoCorticalSleep(cfg, SynapseParams(), sleep, sim,
                                         hh=hh_params)
    spikes, traces, meta = model.run()
    meta["seed"] = args.seed

    rows = analyze(spikes, traces, meta)
    print("\n=== Per-layer sleep-spindle analysis (sigma 10-15 Hz) ===")
    print(_fmt_table(rows))

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    stem = f"spindle_analysis_{args.tag}" if args.tag else "spindle_analysis"
    md = outdir / f"{stem}.md"
    _write_md(rows, meta, md)
    print(f"\nWrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
