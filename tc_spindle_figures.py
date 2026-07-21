"""
Figures demonstrating sleep spindles in the auditory thalamo-cortical model.

Produces a multi-panel figure that shows the spindle from the circuit level up,
in the terms Fernandez & Luthi (2020, Physiol Rev 100:805-868) use:

  (a) a single spindle: LFP + 10-15 Hz band + waxing/waning envelope
  (b) the generating circuit: nRT (RE) bursts alternating with MGB (TC) rebound
      spikes -- the reciprocal loop that *is* the spindle
  (c) several seconds: spindles nested on slow-oscillation UP states
  (d) thalamic spectrogram: discrete spindle events in the sigma band
  (e) infraslow clustering: sigma power oscillating at ~0.02 Hz
      (spindle-rich CONTINUITY vs spindle-poor FRAGILITY periods)
  (f) statistics: intra-spindle frequency and inter-spindle interval

Run (uses the HH model, where spindles are emergent):

    python3 tc_spindle_figures.py --tstop 200000 --outdir out

Long runs are needed for panel (e): >= ~150 s covers 3 infraslow cycles.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from .tc_network import (AuditoryThalamoCorticalSleep, NetworkConfig,
                             SimulationConfig, SynapseParams, SleepParams, HHParams)
    from . import tc_run
except ImportError:
    from tc_network import (AuditoryThalamoCorticalSleep, NetworkConfig,
                            SimulationConfig, SynapseParams, SleepParams, HHParams)
    import tc_run


GREEN = "#2f8f6f"
SIGMA_LO, SIGMA_HI = 10.0, 15.0     # review's spindle band


# ---------------------------------------------------------------------------
#  analysis helpers
# ---------------------------------------------------------------------------

def tc_episodes(mgb_times, gap_ms=500.0):
    """Spindle episodes as contiguous runs of TC (relay) firing.

    This is an ABSOLUTE measure -- unlike an envelope percentile threshold,
    which self-normalises and cannot show spindle suppression.
    Returns (starts, ends, isi_s, dur_s) in ms / s.
    """
    t = np.sort(np.asarray(mgb_times, float))
    if len(t) < 2:
        return np.array([]), np.array([]), np.array([]), np.array([])
    brk = np.where(np.diff(t) > gap_ms)[0]
    starts = np.concatenate(([t[0]], t[brk + 1]))
    ends = np.concatenate((t[brk], [t[-1]]))
    return starts, ends, np.diff(starts) / 1000.0, (ends - starts) / 1000.0


def intra_spindle_freq(mgb_times, nrt_times, starts, ends):
    """Cycles/duration within each episode (review's 'intra-spindle frequency').

    Estimated from the dominant sigma-band peak of the thalamic spike rate
    inside each episode.
    """
    out = []
    allt = np.concatenate([np.asarray(mgb_times, float),
                           np.asarray(nrt_times, float)])
    for s, e in zip(starts, ends):
        if e - s < 400.0:            # too short to estimate a frequency
            continue
        seg = allt[(allt >= s) & (allt <= e)]
        if len(seg) < 20:
            continue
        _, r = tc_run.population_rate(seg - s, e - s, bin_ms=1.0, smooth_ms=3.0)
        pk, _, _, _ = tc_run.detect_peak(r, 1000.0, 7.0, 18.0)
        if np.isfinite(pk):
            out.append(pk)
    return np.array(out)


# ---------------------------------------------------------------------------
#  the figure
# ---------------------------------------------------------------------------

def make_spindle_figure(spikes, traces, meta, out_png, infraslow_hz=0.02):
    tstop = meta["tstop"]
    eeg = tc_run.build_eeg_like(traces, tstop, seed=meta.get("seed", 0))
    t, fs, comp = eeg["t"], eeg["fs"], eeg["eeg"]
    spin = tc_run._bandpass(comp, fs, SIGMA_LO, SIGMA_HI)
    _, env = tc_run.bandpass_envelope(comp, fs, SIGMA_LO, SIGMA_HI)
    slow = tc_run._bandpass(comp, fs, 0.5, 2.0)

    mgb = spikes.get("MGB", {}).get("times", np.array([]))
    nrt = spikes.get("nRT", {}).get("times", np.array([]))
    starts, ends, isi, dur = tc_episodes(mgb)

    fig = plt.figure(figsize=(13, 12))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.0, 1.0, 0.9], hspace=0.55,
                          wspace=0.22)

    # ---- pick a representative spindle for the zoom panels ----
    if len(starts):
        k = int(np.argmax(dur[:len(starts)])) if len(dur) else 0
        z0, z1 = starts[k] - 400.0, min(starts[k] + 2200.0, tstop)
    else:
        z0, z1 = 0.0, min(3000.0, tstop)
    mz = (t >= z0) & (t <= z1)

    # ---- (a) single spindle: LFP, sigma band, envelope ----
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(t[mz], comp[mz], color="0.55", lw=0.8, label="LFP (composite)")
    ax.plot(t[mz], spin[mz], color=GREEN, lw=1.1,
            label=f"spindle band {SIGMA_LO:.0f}-{SIGMA_HI:.0f} Hz")
    ax.plot(t[mz], env[mz], color="k", lw=1.0, alpha=0.65, label="envelope")
    ax.plot(t[mz], -env[mz], color="k", lw=1.0, alpha=0.65)
    ax.set_title("(a) A single spindle: waxing and waning", fontsize=10, loc="left")
    ax.set_xlabel("time (ms)"); ax.set_ylabel("amplitude (uV)")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    # ---- (b) the generating circuit: RE bursts vs TC rebound ----
    ax = fig.add_subplot(gs[0, 1])
    for times, y, c, lab in [(nrt, 1.0, "#2471a3", "nRT / RE (inhibitory)"),
                             (mgb, 0.0, "#c0392b", "MGB / TC (relay)")]:
        sel = times[(times >= z0) & (times <= z1)]
        ax.plot(sel, np.full_like(sel, y), "|", ms=11, color=c, label=lab)
    ax.set_ylim(-0.6, 1.6); ax.set_yticks([0, 1])
    ax.set_yticklabels(["TC", "RE"])
    ax.set_title("(b) The RE<->TC loop that generates it", fontsize=10, loc="left")
    ax.set_xlabel("time (ms)")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    # ---- (c) spindles nested on slow-oscillation UP states ----
    ax = fig.add_subplot(gs[1, :])
    w1 = min(z0 + 12000.0, tstop)
    mw = (t >= z0) & (t <= w1)
    ax.plot(t[mw], slow[mw], color="k", lw=1.6, label="slow wave 0.5-2 Hz")
    ax.plot(t[mw], spin[mw] * 1.0, color=GREEN, lw=0.8, label="spindle band")
    for s, e in zip(starts, ends):
        if e > z0 and s < w1:
            ax.axvspan(s, e, color=GREEN, alpha=0.16, lw=0)
    ax.set_title("(c) Spindles (shaded) nested on slow-oscillation UP states",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (ms)"); ax.set_ylabel("amplitude (uV)")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    # ---- (d) thalamic spectrogram: discrete spindle events ----
    ax = fig.add_subplot(gs[2, :])
    allt = np.concatenate([mgb, nrt]) if len(mgb) + len(nrt) else np.array([0.0])
    _, rate = tc_run.population_rate(allt, tstop, bin_ms=1.0, smooth_ms=3.0)
    seg = min(2048, max(256, len(rate) // 40))
    Pxx, freqs_s, bins_s, im = ax.specgram(
        rate - rate.mean(), NFFT=seg, Fs=1000.0,
        noverlap=int(seg * 0.85), cmap="magma")
    # Contrast-stretch to the 60th-99.5th percentile of in-band power, otherwise
    # the broadband 1/f background washes the spindle events out.
    band = (freqs_s >= 2.0) & (freqs_s <= 25.0)
    dB = 10.0 * np.log10(np.maximum(Pxx[band], 1e-20))
    im.set_clim(np.percentile(dB, 60), np.percentile(dB, 99.5))
    ax.axhline(SIGMA_LO, color="w", ls="--", lw=0.8, alpha=0.7)
    ax.axhline(SIGMA_HI, color="w", ls="--", lw=0.8, alpha=0.7)
    ax.set_ylim(0, 25)
    ax.set_title("(d) Thalamic spectrogram: discrete events in the sigma band "
                 f"({SIGMA_LO:.0f}-{SIGMA_HI:.0f} Hz, dashed)",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (s)"); ax.set_ylabel("frequency (Hz)")

    # ---- (e) infraslow clustering of sigma power ----
    ax = fig.add_subplot(gs[3, 0])
    need_s = 3.0 / max(1e-6, infraslow_hz)
    if tstop / 1000.0 >= need_s:
        # smooth the envelope to expose the infraslow (~50 s) rhythm
        win = max(3, int(5.0 * fs))                     # ~5 s boxcar
        sm = np.convolve(env, np.ones(win) / win, mode="same")
        ax.plot(t / 1000.0, sm, color=GREEN, lw=1.0)
        e2 = sm[int(2 * fs):] - sm[int(2 * fs):].mean()
        f = np.fft.rfftfreq(len(e2), 1.0 / fs)
        P = np.abs(np.fft.rfft(e2)) ** 2
        band = (f > infraslow_hz * 0.6) & (f < infraslow_hz * 1.4)
        ref = (f > 0.01) & (f < 1.0)
        if band.any() and ref.any():
            pk = f[band][np.argmax(P[band])]
            ratio = P[band].max() / max(1e-12, np.median(P[ref]))
            ax.set_title(f"(e) Infraslow clustering: sigma power at "
                         f"{pk:.3f} Hz (ratio {ratio:.1f})",
                         fontsize=10, loc="left")
        ax.set_xlabel("time (s)"); ax.set_ylabel("sigma envelope (uV)")
    else:
        ax.text(0.5, 0.5, f"(e) needs a run >= {need_s:.0f} s\n"
                          f"to resolve {infraslow_hz} Hz clustering",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_axis_off()

    # ---- (f) statistics ----
    ax = fig.add_subplot(gs[3, 1])
    freqs = intra_spindle_freq(mgb, nrt, starts, ends)
    txt = []
    if len(freqs):
        ax.hist(freqs, bins=np.arange(6, 19, 1.0), color=GREEN, alpha=0.85,
                edgecolor="k", lw=0.5)
        txt.append(f"intra-spindle freq: {np.median(freqs):.1f} Hz")
    ax.axvspan(SIGMA_LO, SIGMA_HI, color="k", alpha=0.08, lw=0)
    if len(isi):
        txt.append(f"inter-spindle interval: {np.median(isi):.1f} s")
    if len(dur):
        txt.append(f"duration: {np.median(dur):.2f} s")
    if len(starts):
        txt.append(f"density: {len(starts) / (tstop / 60000.0):.1f} /min")
    ax.set_title("(f) Intra-spindle frequency (shaded = review 10-15 Hz)",
                 fontsize=10, loc="left")
    ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("count")
    if txt:
        ax.text(0.98, 0.97, "\n".join(txt), transform=ax.transAxes, fontsize=8,
                va="top", ha="right",
                bbox=dict(boxstyle="round", fc="w", ec="0.7", alpha=0.9))

    fig.suptitle("Sleep spindles in the auditory thalamo-cortical model "
                 "(ht_neuron, emergent)", fontsize=13, y=0.98)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved spindle figure to {out_png}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=str,
                    default="config/network_auditory_local.yaml")
    ap.add_argument("--tstop", type=float, default=200000.0,
                    help="ms; >= 150000 to resolve the 0.02 Hz clustering")
    ap.add_argument("--outdir", type=str, default="out")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-trigger", action="store_true",
                    help="disable the external modulatory spindle trigger")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    cfg = NetworkConfig.from_file(args.config)
    cfg.tstop = args.tstop
    sim = SimulationConfig(seed=args.seed, neuron_model="ht_neuron",
                           record_traces=True)
    sleep = SleepParams(emergent_spindles=True,
                        spindle_trigger=not args.no_trigger)
    print(f"Running HH model for spindle figures ({cfg.tstop/1000:.0f} s)...")
    model = AuditoryThalamoCorticalSleep(cfg, SynapseParams(), sleep, sim)
    spikes, traces, meta = model.run()
    meta["seed"] = args.seed
    make_spindle_figure(spikes, traces, meta,
                        outdir / "tc_spindles.png",
                        infraslow_hz=HHParams().infraslow_freq)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
