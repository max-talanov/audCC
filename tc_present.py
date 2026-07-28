"""
Presentation-quality "indicative picture" of slow waves and sleep spindles,
for demonstrating the auditory thalamo-cortical model to a neuroscientist.

Renders one clean figure in the canonical sleep-EEG style (cf. Fernandez &
Luthi 2020; the Reuniens LFP / SP-SW decomposition figures): a composite
auditory-cortex LFP with sleep spindles nested on slow-oscillation UP states,
decomposed beneath into its Slow oscillation (0.5-2 Hz) and Sleep spindle
(10-15 Hz) bands, with scale bars and annotations.

    python3 tc_present.py --config config/network_auditory_mushtaq.yaml \
        --tstop 12000 --outdir out
"""

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    from tc_sleep.tc_network import (AuditoryThalamoCorticalSleep, NetworkConfig,
                                     SimulationConfig, SynapseParams, SleepParams)
    from tc_sleep import tc_run
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tc_network import (AuditoryThalamoCorticalSleep, NetworkConfig,
                            SimulationConfig, SynapseParams, SleepParams)
    import tc_run

LFP = "#334155"      # slate
SPIN = "#0d9488"     # teal
SLOW = "#1e3a8a"     # navy


def _scalebar(ax, x0, y0, dx, dy, xlab, ylab):
    ax.plot([x0, x0], [y0, y0 + dy], color="k", lw=1.6)
    ax.plot([x0, x0 + dx], [y0, y0], color="k", lw=1.6)
    xr = ax.get_xlim()[1] - ax.get_xlim()[0]
    yr = ax.get_ylim()[1] - ax.get_ylim()[0]
    ax.text(x0 - 0.008 * xr, y0 + dy / 2, ylab, ha="right", va="center", fontsize=9)
    ax.text(x0 + dx / 2, y0 - 0.05 * yr, xlab, ha="center", va="top", fontsize=9)


def make_figure(spikes, traces, meta, out_png, window_s=5.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tstop = meta["tstop"]
    eeg = tc_run.build_eeg_like(traces, tstop, seed=meta.get("seed", 0))
    t, fs, comp = eeg["t"], eeg["fs"], eeg["eeg"]
    spin = tc_run._bandpass(comp, fs, 10.0, 15.0)
    _, env = tc_run.bandpass_envelope(comp, fs, 10.0, 15.0)
    sp_win, _ = tc_run.detect_sp_sw(eeg)

    # Slow oscillation derived from CORTICAL FIRING RATE (unambiguous UP/DOWN:
    # UP = cortex active), using whichever cortical layers are active.
    cort = [l for l in ["L23", "L5", "L6", "L4"] if l in spikes
            and len(spikes[l]["times"]) > 20]
    crate = None
    if cort:
        allc = np.concatenate([spikes[l]["times"] for l in cort])
        _, crate = tc_run.population_rate(allc, tstop, bin_ms=1.0, smooth_ms=25.0)
        slow = tc_run._bandpass(crate, fs, 0.5, 2.0)
    else:
        slow = tc_run._bandpass(comp, fs, 0.5, 2.0)

    # Measure the spindle-to-SO coupling from the data and LABEL the figure to
    # match (never assert "UP-nested" unless the cortex is actually active at
    # the spindle). ratio = cortical firing at spindle peaks / overall mean.
    from scipy.signal import find_peaks
    coupling_label = "phase-locked to the slow oscillation (~1 per cycle)"
    ratio = np.nan
    if crate is not None:
        cr = crate[:len(env)]
        pk, _ = find_peaks(env, height=np.percentile(env, 80), distance=300)
        if len(pk) and cr.mean() > 0:
            ratio = float(cr[pk].mean() / cr.mean())
            if ratio > 1.1:
                coupling_label = "nested on the slow-oscillation UP state"
            elif ratio < 0.9:
                coupling_label = ("phase-locked to the slow oscillation "
                                  "(in the DOWN phase)")
    print(f"Spindle-SO coupling: cortical firing at spindle / mean = "
          f"{ratio:.2f}  -> {coupling_label}")

    # choose a clean window: start at a detected spindle, span window_s
    w0 = max(500.0, (sp_win[0][0] - 900.0) if sp_win else 500.0)
    w1 = min(w0 + window_s * 1000.0, tstop)
    m = (t >= w0) & (t <= w1)
    tt = (t[m] - w0) / 1000.0                      # seconds, from 0

    fig, axes = plt.subplots(3, 1, figsize=(11, 6.6), sharex=True)
    fig.subplots_adjust(hspace=0.4, left=0.06, right=0.97, top=0.9, bottom=0.08)

    def clean(ax):
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])

    # slow wave scaled to the composite LFP range, for the overlay
    sw = slow[m]
    sw = sw / (np.percentile(np.abs(sw), 98) or 1.0) * (0.55 * comp[m].std() * 2)

    # ---- (top) composite LFP + slow-wave overlay, spindles shaded ----
    ax = axes[0]
    for (a, b) in sp_win:
        if b > w0 and a < w1:
            ax.axvspan((a - w0) / 1000.0, (b - w0) / 1000.0,
                       color=SPIN, alpha=0.10, lw=0)
    ax.plot(tt, comp[m], color=LFP, lw=0.8, label="LFP (composite)")
    ax.plot(tt, sw, color=SLOW, lw=2.0, alpha=0.85, label="slow oscillation")
    ax.set_title(f"Auditory-cortex LFP  —  sleep spindles (shaded) "
                 f"{coupling_label}", fontsize=11, loc="left")
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9, ncol=2)
    clean(ax)
    yr = comp[m].max() - comp[m].min()
    _scalebar(ax, tt[-1] - 0.55, comp[m].min(), 0.5, 0.5 * yr,
              "0.5 s", f"{0.5 * yr:.0f} µV")

    # ---- (middle) spindle band ----
    ax = axes[1]
    ax.plot(tt, spin[m], color=SPIN, lw=0.9)
    ax.plot(tt, env[m], color="k", lw=0.8, alpha=0.5)
    ax.plot(tt, -env[m], color="k", lw=0.8, alpha=0.5)
    ax.set_title("Sleep spindle band (10–15 Hz)  —  discrete, "
                 "waxing/waning", fontsize=11, loc="left", color=SPIN)
    clean(ax)
    sr = spin[m].max() - spin[m].min()
    _scalebar(ax, tt[-1] - 0.55, spin[m].min(), 0.5, 0.6 * sr,
              "0.5 s", f"{0.6 * sr:.0f} µV")

    # ---- (bottom) slow oscillation (normalised; it is a firing-rate proxy) ----
    ax = axes[2]
    swn = slow[m] / (np.percentile(np.abs(slow[m]), 98) or 1.0)
    ax.plot(tt, swn, color=SLOW, lw=1.6)
    ax.set_title("Slow oscillation (0.5–2 Hz, from cortical firing)  —  "
                 "UP (peak) / DOWN (trough) states",
                 fontsize=11, loc="left", color=SLOW)
    clean(ax)
    _scalebar(ax, tt[-1] - 0.55, swn.min(), 0.5, 1.0, "0.5 s", "cortical\nfiring (a.u.)")

    fig.suptitle("Slow oscillations and sleep spindles — auditory "
                 "thalamo-cortical model (NEST)", fontsize=13, y=0.975)
    fig.savefig(out_png, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved presentation figure to {out_png}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=str,
                    default="config/network_auditory_adex.yaml")
    ap.add_argument("--tstop", type=float, default=12000.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", type=str, default="out")
    ap.add_argument("--window", type=float, default=5.0, help="window shown (s)")
    args = ap.parse_args(argv)

    cfg = NetworkConfig.from_file(args.config); cfg.tstop = args.tstop
    nm = tc_run._config_neuron_model(args.config) or "iaf_cond_exp"
    sim = SimulationConfig(seed=args.seed, neuron_model=nm, record_traces=True)
    sleep = (SleepParams(emergent_spindles=True) if nm != "iaf_cond_exp"
             else SleepParams())
    print(f"Running column ({args.tstop/1000:.0f} s, model {nm})...")
    model = AuditoryThalamoCorticalSleep(cfg, SynapseParams(), sleep, sim)
    spikes, traces, meta = model.run(); meta["seed"] = args.seed

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    make_figure(spikes, traces, meta,
                outdir / "slow_waves_and_spindles.png", window_s=args.window)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
