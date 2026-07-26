"""
Figures illustrating sleep spindles in the NEURON thalamo-cortical model.

Shows the spindle from the single-cell mechanism up to the network oscillation,
in the terms of Fernandez & Luthi (2020):

  (a) TC relay cell: a low-threshold Ca2+ rebound BURST (Destexhe I_T)
  (b) RE reticular cell: rebound burst (Ca_v3.3 I_T + SK2)
  (c) gap-junction synchronisation of two RE cells (coupled vs not)
  (d) network spike raster (TC + RE), showing the reciprocal loop
  (e) thalamic LFP-proxy + 8-15 Hz spindle-band filter
  (f) spectrogram of the thalamic signal (the ~10 Hz spindle-band oscillation)

Run (from the neuron/ directory, in the migration venv):

    ../.venv-neuron/bin/python tc_neuron_figures.py --outdir ../out
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, hilbert
from neuron import h

try:
    from . import tc_neuron as T
    from . import tc_network_nrn as N
except ImportError:
    import tc_neuron as T
    import tc_network_nrn as N

GREEN = "#2f8f6f"
RED = "#c0392b"
BLUE = "#2471a3"


def _record_vm(cell):
    return np.asarray(cell.t), np.asarray(cell.vsoma)


def single_cell_rebound(make_cell, amp=-0.2, dur=500.0):
    """Hyperpolarise, release, return (t, Vm, spike count) of the rebound burst."""
    cell = make_cell()
    cell.record()
    ic = h.IClamp(cell.soma(0.5)); ic.delay, ic.dur, ic.amp = 300, dur, amp
    h.celsius = 36
    h.finitialize(-74 if isinstance(cell, T.TCCell) else -75)
    h.continuerun(300 + dur + 250)
    t, v = _record_vm(cell)
    rel = 300 + dur
    sp = np.asarray(cell.spikes)
    n = int(((sp > rel) & (sp < rel + 250)).sum())
    return t, v, n, rel


def gap_pair(g_gap):
    """Drive RE cell A only; return (t, Vm_A, Vm_B) to show electrical coupling."""
    A, B = T.RECell(gsk=0.0), T.RECell(gsk=0.0)
    A.record(); B.record()
    if g_gap > 0:
        T.gap_junction(A, B, g=g_gap)
    h.celsius = 36
    h.finitialize(-64); h.continuerun(150)
    ic = h.IClamp(A.soma(0.5)); ic.delay, ic.dur, ic.amp = 150, 200, 0.25
    h.continuerun(450)
    ta, va = _record_vm(A); tb, vb = _record_vm(B)
    return ta, va, vb


def _bandpass(x, fs, lo, hi):
    x = np.asarray(x, float) - np.mean(x)
    b, a = butter(3, [lo / (fs / 2), hi / (fs / 2)], btype="band")
    return filtfilt(b, a, x)


def make_figure(out_png, net_tstop=8000.0):
    fig = plt.figure(figsize=(13, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.1], hspace=0.42,
                          wspace=0.22)

    # ---- (a) TC rebound burst ----
    ax = fig.add_subplot(gs[0, 0])
    t, v, n, rel = single_cell_rebound(lambda: T.TCCell(gsk=0.0))
    m = (t > rel - 120) & (t < rel + 160)
    ax.plot(t[m] - rel, v[m], color=GREEN, lw=0.9)
    ax.axvline(0, color="0.6", ls="--", lw=0.8)
    ax.set_title(f"(a) TC relay cell: {n}-spike rebound burst (Destexhe I_T)",
                 fontsize=10, loc="left")
    ax.set_xlabel("time from release (ms)"); ax.set_ylabel("V_m (mV)")

    # ---- (b) RE rebound burst ----
    ax = fig.add_subplot(gs[0, 1])
    t, v, n, rel = single_cell_rebound(lambda: T.RECell(), amp=-0.2)
    m = (t > rel - 120) & (t < rel + 160)
    ax.plot(t[m] - rel, v[m], color=BLUE, lw=0.9)
    ax.axvline(0, color="0.6", ls="--", lw=0.8)
    ax.set_title(f"(b) RE reticular cell: {n}-spike burst (Ca_v3.3 I_T + SK2)",
                 fontsize=10, loc="left")
    ax.set_xlabel("time from release (ms)"); ax.set_ylabel("V_m (mV)")

    # ---- (c) gap-junction synchronisation ----
    ax = fig.add_subplot(gs[1, 0])
    ta, va, vb0 = gap_pair(0.0)      # uncoupled B
    _, _, vb1 = gap_pair(0.02)       # coupled B
    m = (ta > 140) & (ta < 360)
    ax.plot(ta[m], va[m], color="0.4", lw=0.8, label="cell A (driven)")
    ax.plot(ta[m], vb0[m] + 130, color="0.7", lw=0.8, label="cell B, g=0 (silent)")
    ax.plot(ta[m], vb1[m] + 260, color=RED, lw=0.8, label="cell B, g=0.02 (recruited)")
    ax.set_title("(c) Gap-junction synchronisation of two RE cells",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (ms)"); ax.set_yticks([])
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    # ---- run the network for the remaining panels ----
    net = N.ThalamicNet()
    spikes, traces, meta = net.run(tstop=net_tstop)
    grid = traces["MGB"]["time"]
    thal_v = traces["MGB"]["voltage"]

    # ---- (d) raster ----
    ax = fig.add_subplot(gs[1, 1])
    tc = spikes["MGB"]; re = spikes["nRT"]
    ax.plot(re["times"], re["senders"], "|", ms=5, color=BLUE, label="nRT (RE)")
    ax.plot(tc["times"], tc["senders"], "|", ms=5, color=RED, label="MGB (TC)")
    ax.set_title("(d) Network spike raster (RE inhibits TC, TC excites RE)",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (ms)"); ax.set_ylabel("cell #")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    # ---- (e) thalamic LFP proxy + spindle band ----
    ax = fig.add_subplot(gs[2, 0])
    fs = 1000.0
    lfp = thal_v - thal_v.mean()
    spin = _bandpass(lfp, fs, 8.0, 15.0)
    env = np.abs(hilbert(spin))
    w = (grid > 500) & (grid < min(4500, net_tstop))
    ax.plot(grid[w], lfp[w], color="0.6", lw=0.6, label="thalamic LFP proxy")
    ax.plot(grid[w], spin[w], color=GREEN, lw=0.9, label="8-15 Hz spindle band")
    ax.plot(grid[w], env[w], color="k", lw=0.8, alpha=0.6)
    ax.plot(grid[w], -env[w], color="k", lw=0.8, alpha=0.6)
    ax.set_title("(e) Thalamic LFP proxy and spindle-band (8-15 Hz) component",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (ms)"); ax.set_ylabel("mV")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    # ---- (f) spectrogram ----
    ax = fig.add_subplot(gs[2, 1])
    allt = np.concatenate([tc["times"], re["times"]])
    rate = np.histogram(allt, bins=np.arange(0, net_tstop + 1, 1.0))[0].astype(float)
    seg = 1024
    Pxx, fr, bins, im = ax.specgram(rate - rate.mean(), NFFT=seg, Fs=1000.0,
                                    noverlap=int(seg * 0.9), cmap="magma")
    band = (fr >= 2) & (fr <= 25)
    dB = 10 * np.log10(np.maximum(Pxx[band], 1e-20))
    im.set_clim(np.percentile(dB, 60), np.percentile(dB, 99.5))
    ax.axhline(8, color="w", ls="--", lw=0.7, alpha=0.7)
    ax.axhline(15, color="w", ls="--", lw=0.7, alpha=0.7)
    ax.set_ylim(0, 25)
    ax.set_title("(f) Thalamic spectrogram (spindle band 8-15 Hz dashed)",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (s)"); ax.set_ylabel("frequency (Hz)")

    fig.suptitle("Sleep spindles in the NEURON thalamo-cortical model "
                 "(conductance-based, mechanistic)", fontsize=13, y=0.985)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved spindle figure to {out_png}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=str, default="../out")
    ap.add_argument("--tstop", type=float, default=8000.0)
    args = ap.parse_args(argv)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    make_figure(outdir / "tc_neuron_spindles.png", net_tstop=args.tstop)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
