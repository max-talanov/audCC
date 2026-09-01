"""
Post-hoc analysis figures for a ctx_thalamus_mpi.py run (loads the saved
.npz spike dump -- does not need NEURON or MPI).

Two figures:

  1. <tag>_meanfield.png -- population-rate mean-field proxy (thalamus vs.
     cortex), slow-oscillation band (0.5-4 Hz), spindle band (8-15 Hz) with
     Hilbert envelope, and a spectrogram of the thalamic (TC+RE) signal.

  2. <tag>_spindles.png -- RE burst shape (raster, ISI histogram, burst-size
     histogram) and a spindle-event-triggered PSTH across TC/RE/L4E/L23E/L5E/L6E
     showing whether/how a thalamic spindle volley propagates up the column.

Usage:
    python3 neuron/ctx_analyze.py res/2026-08-31/ctx_nrn_45171023.npz --outdir out
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, hilbert

GREEN = "#2f8f6f"
RED = "#c0392b"
BLUE = "#2471a3"
GREY = "#7f8c8d"


def _bandpass(x, fs, lo, hi, order=3):
    x = np.asarray(x, float) - np.mean(x)
    sos = butter(order, [lo / (fs / 2), hi / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x)


def _rate(times, tstop, bin_ms=1.0, smooth_ms=0.0):
    bins = np.arange(0, tstop + bin_ms, bin_ms)
    r = np.histogram(times, bins=bins)[0].astype(float)
    if smooth_ms > 0:
        k = max(1, int(round(smooth_ms / bin_ms)))
        r = np.convolve(r, np.ones(k) / k, mode="same")
    return r, bins[:-1]


def _re_burst_stats(t, g, re_lo, re_hi, burst_gap=30.0, event_gap=300.0):
    """Same definition as ctx_thalamus_mpi.py's _re_burst_stats (kept in sync
    manually since that module cannot be imported without NEURON/MPI on the
    plotting host), plus the raw ISI/burst-size arrays for histograms."""
    m = (g >= re_lo) & (g < re_hi)
    tp, gp = t[m], g[m]
    if len(tp) < 2:
        return dict(n_spikes=len(tp), frac_burst=0.0, mean_burst_size=0.0,
                    n_events=0, event_hz=0.0, isis=np.array([]),
                    burst_sizes=np.array([]))
    order = np.argsort(gp, kind="stable")
    tp, gp = tp[order], gp[order]
    same_cell = np.diff(gp) == 0
    isi = np.diff(tp)
    in_burst = same_cell & (isi < burst_gap)
    frac_burst = float(in_burst.sum()) / len(tp)

    burst_sizes = []
    run = 1
    for i in range(len(in_burst)):
        if in_burst[i]:
            run += 1
        else:
            if run > 1:
                burst_sizes.append(run)
            run = 1
    if run > 1:
        burst_sizes.append(run)
    mean_burst = float(np.mean(burst_sizes)) if burst_sizes else 0.0

    tsort = np.sort(t[m])
    starts = tsort[np.diff(tsort, prepend=-1e9) > event_gap]
    event_hz = 1000.0 / np.mean(np.diff(starts)) if len(starts) > 1 else 0.0

    isis_within_cell = isi[same_cell]

    return dict(n_spikes=len(tp), frac_burst=round(frac_burst, 3),
                mean_burst_size=round(mean_burst, 2), n_events=len(starts),
                event_hz=round(event_hz, 3), isis=isis_within_cell,
                burst_sizes=np.array(burst_sizes))


def make_meanfield_figure(npz, out_png, window=(20000, 30000)):
    t, g, ranges, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    fs = 1000.0  # 1 ms bins

    tc_lo, tc_hi = ranges["tc"]
    re_lo, re_hi = ranges["re"]
    thal = t[(g >= tc_lo) & (g < re_hi)]
    cx_layers = ["l4e", "l23e", "l5e", "l6e"]
    cx = t[np.isin(g, np.concatenate([np.arange(*ranges[k]) for k in cx_layers]))]

    thal_rate, tbins = _rate(thal, tstop, bin_ms=1.0, smooth_ms=5.0)
    cx_rate, _ = _rate(cx, tstop, bin_ms=1.0, smooth_ms=5.0)

    thal_so = _bandpass(thal_rate, fs, 0.5, 4.0)
    thal_sp = _bandpass(thal_rate, fs, 8.0, 15.0)
    thal_env = np.abs(hilbert(thal_sp))
    cx_so = _bandpass(cx_rate, fs, 0.5, 4.0)

    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2], hspace=0.45, wspace=0.25)

    ax = fig.add_subplot(gs[0, :])
    ax.plot(tbins / 1000.0, thal_rate, color=BLUE, lw=0.5, label="thalamus (TC+RE) rate")
    ax.plot(tbins / 1000.0, cx_rate, color=RED, lw=0.5, alpha=0.7, label="cortex (L4/L2-3/L5/L6 E) rate")
    ax.set_title("(a) Population firing-rate mean-field, full 200 s run", fontsize=10, loc="left")
    ax.set_xlabel("time (s)"); ax.set_ylabel("spikes / ms (1 ms bin, 5 ms smoothed)")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    w = (tbins >= window[0]) & (tbins < window[1])
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(tbins[w] / 1000.0, thal_so[w], color=GREEN, lw=1.0, label="thalamic 0.5-4 Hz (SO)")
    ax.plot(tbins[w] / 1000.0, cx_so[w], color=GREY, lw=1.0, label="cortical 0.5-4 Hz (SO)")
    ax.set_title(f"(b) Slow-oscillation band, {window[0]/1000:.0f}-{window[1]/1000:.0f} s window",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (s)"); ax.set_ylabel("a.u.")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    ax = fig.add_subplot(gs[1, 1])
    ax.plot(tbins[w] / 1000.0, thal_sp[w], color=BLUE, lw=0.9, label="8-15 Hz spindle band")
    ax.plot(tbins[w] / 1000.0, thal_env[w], color="k", lw=0.8, alpha=0.6)
    ax.plot(tbins[w] / 1000.0, -thal_env[w], color="k", lw=0.8, alpha=0.6)
    ax.set_title("(c) Thalamic spindle-band component + envelope", fontsize=10, loc="left")
    ax.set_xlabel("time (s)"); ax.set_ylabel("a.u.")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    ax = fig.add_subplot(gs[2, :])
    seg = 2048
    Pxx, fr, bins_, im = ax.specgram(thal_rate - thal_rate.mean(), NFFT=seg, Fs=fs,
                                     noverlap=int(seg * 0.9), cmap="magma")
    band = (fr >= 0.2) & (fr <= 20)
    dB = 10 * np.log10(np.maximum(Pxx[band], 1e-20))
    im.set_clim(np.percentile(dB, 55), np.percentile(dB, 99.5))
    ax.axhline(0.5, color="c", ls="--", lw=0.7, alpha=0.7)
    ax.axhline(4, color="c", ls="--", lw=0.7, alpha=0.7)
    ax.axhline(8, color="w", ls="--", lw=0.7, alpha=0.7)
    ax.axhline(15, color="w", ls="--", lw=0.7, alpha=0.7)
    ax.set_ylim(0, 20)
    ax.set_title("(d) Thalamic spectrogram (SO band cyan dashed, spindle band white dashed)",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (s)"); ax.set_ylabel("frequency (Hz)")

    fig.suptitle("Mean-field slow-oscillation / spindle indication -- full corticothalamic network",
                 fontsize=13, y=0.99)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved mean-field figure to {out_png}")


def make_spindle_figure(npz, out_png, window=(20000, 30000)):
    t, g, ranges, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    re_lo, re_hi = ranges["re"]
    tc_lo, tc_hi = ranges["tc"]

    s = _re_burst_stats(t, g, re_lo, re_hi)

    fig = plt.figure(figsize=(13, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.1], hspace=0.48, wspace=0.25)

    w = (t >= window[0]) & (t < window[1])
    ax = fig.add_subplot(gs[0, :])
    for name, color in [("tc", BLUE), ("re", RED)]:
        lo, hi = ranges[name]
        m = w & (g >= lo) & (g < hi)
        ax.plot(t[m] / 1000.0, g[m], "|", ms=3, color=color, label=name.upper())
    ax.set_title(f"(a) TC/RE raster, {window[0]/1000:.0f}-{window[1]/1000:.0f} s "
                 f"(RE burst frac={s['frac_burst']:.2f}, mean burst size={s['mean_burst_size']:.1f})",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (s)"); ax.set_ylabel("gid")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    ax = fig.add_subplot(gs[1, 0])
    if len(s["isis"]):
        ax.hist(s["isis"], bins=np.arange(0, 205, 5), color=BLUE, alpha=0.85)
    ax.axvline(30, color=RED, ls="--", lw=1.0, label="burst-ISI threshold (30 ms)")
    ax.set_title("(b) RE within-cell ISI distribution", fontsize=10, loc="left")
    ax.set_xlabel("ISI (ms)"); ax.set_ylabel("count")
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    bs = s["burst_sizes"]
    if len(bs):
        ax.hist(bs, bins=np.arange(0.5, bs.max() + 1.5, 1), color=GREEN, alpha=0.85)
    ax.set_title(f"(c) RE burst-size distribution (n_events={s['n_events']}, "
                 f"event_Hz={s['event_hz']:.2f})", fontsize=10, loc="left")
    ax.set_xlabel("spikes per burst"); ax.set_ylabel("count")

    # (d) event-triggered PSTH: align on each RE population volley, show
    # per-layer rate around it to see thalamus->cortex propagation.
    t_re = np.sort(t[(g >= re_lo) & (g < re_hi)])
    if len(t_re) > 1:
        gaps = np.diff(t_re)
        ev_times = t_re[np.concatenate([[True], gaps > 300.0])]
    else:
        ev_times = np.array([])

    layers = ["tc", "re", "l4e", "l23e", "l5e", "l6e"]
    colors = [BLUE, RED, "#8e44ad", "#d35400", GREEN, GREY]
    pre, post, bin_ms = 200.0, 400.0, 5.0
    edges = np.arange(-pre, post + bin_ms, bin_ms)
    ax = fig.add_subplot(gs[2, :])
    n_ev = min(len(ev_times), 200)
    ev_sub = ev_times[np.linspace(0, len(ev_times) - 1, n_ev).astype(int)] if len(ev_times) else ev_times
    for name, color in zip(layers, colors):
        lo, hi = ranges[name]
        m = (g >= lo) & (g < hi)
        tl, gl = t[m], g[m]
        n_cells = hi - lo
        counts = np.zeros(len(edges) - 1)
        for ev in ev_sub:
            rel = tl - ev
            keep = (rel >= -pre) & (rel < post)
            counts += np.histogram(rel[keep], bins=edges)[0]
        if n_ev > 0 and n_cells > 0:
            rate_hz = counts / n_ev / n_cells / (bin_ms / 1000.0)
        else:
            rate_hz = counts
        ax.plot(edges[:-1] + bin_ms / 2, rate_hz, color=color, lw=1.2, label=name.upper())
    ax.axvline(0, color="k", ls="--", lw=0.8, alpha=0.6)
    ax.set_title(f"(d) RE-volley-triggered per-layer PSTH (n={n_ev} events) -- "
                 "thalamus -> cortex propagation", fontsize=10, loc="left")
    ax.set_xlabel("time from RE volley onset (ms)"); ax.set_ylabel("rate (Hz/cell)")
    ax.legend(fontsize=8, ncol=3, loc="upper right", framealpha=0.9)

    fig.suptitle("Spindle / RE-burst-shape analysis -- full corticothalamic network",
                 fontsize=13, y=0.99)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved spindle-analysis figure to {out_png}")
    return s


def _synaptic_lfp(times, tstop, fs=1000.0, tau_rise=1.0, tau_decay=10.0, kernel_ms=60.0):
    """Convolve a population spike train with an alpha-like PSP kernel to get
    a smooth current-like proxy (closer to a real LFP than raw binned rate --
    a real LFP is dominated by synaptic currents, which are much slower and
    smoother than the spikes that trigger them)."""
    bin_ms = 1000.0 / fs
    bins = np.arange(0, tstop + bin_ms, bin_ms)
    counts = np.histogram(times, bins=bins)[0].astype(float)
    tk = np.arange(0, kernel_ms, bin_ms)
    kernel = (np.exp(-tk / tau_decay) - np.exp(-tk / tau_rise))
    kernel /= kernel.max() if kernel.max() > 0 else 1.0
    lfp = np.convolve(counts, kernel, mode="full")[:len(counts)]
    return lfp, bins[:-1]


def make_lfp_figure(npz, out_png, window=(20000, 30000), epoch_ms=2000.0):
    """LFP-style figure for direct comparison to spindle-review figures
    (e.g. Fernandez & Luthi 2020): stacked laminar traces, a single-epoch
    zoom, and a composite-LFP raw/spindle-band/spectrogram panel."""
    t, g, ranges, tstop = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    fs = 1000.0

    # Laminar order, superficial to deep, cortex only (thalamus shown separately)
    lam_layers = ["l23e", "l4e", "l5e", "l6e"]
    lam_labels = ["L2/3", "L4", "L5", "L6"]
    lam_colors = ["#1b4965", "#2471a3", "#5fa8d3", "#7f8c8d"]

    lfps = {}
    for name in lam_layers + ["tc", "re"]:
        lo, hi = ranges[name]
        m = (g >= lo) & (g < hi)
        lfp, bins = _synaptic_lfp(t[m], tstop, fs=fs)
        # sign convention: excitatory synaptic current -> deflection; deep
        # layers flipped to mimic the classic surface-negative/deep-positive
        # dipole seen in real laminar LFP during synchronized volleys
        sign = -1.0 if name in ("l5e", "l6e", "re") else 1.0
        lfps[name] = sign * (lfp - lfp.mean()) / (lfp.std() + 1e-9)

    composite = np.mean([lfps[k] for k in lam_layers], axis=0)
    thal = np.mean([lfps["tc"], lfps["re"]], axis=0)

    fig = plt.figure(figsize=(13, 12))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.3, 1, 1.2], hspace=0.5)

    # ---- (a) stacked laminar + thalamic traces, literature-style ----
    ax = fig.add_subplot(gs[0])
    w = (bins >= window[0]) & (bins < window[1])
    tsec = bins[w] / 1000.0
    offset = 0.0
    step = 5.0
    for name, label, color in zip(["re", "tc"] + lam_layers,
                                   ["RE (thalamic)", "TC (thalamic)"] + lam_labels,
                                   ["#c0392b", "#2f8f6f"] + lam_colors):
        ax.plot(tsec, lfps[name][w] + offset, color=color, lw=0.7)
        ax.text(tsec[0] - 0.15, offset, label, ha="right", va="center", fontsize=8)
        offset -= step
    ax.set_yticks([])
    ax.set_xlabel("time (s)")
    ax.set_title(f"(a) Laminar LFP proxy, {window[0]/1000:.0f}-{window[1]/1000:.0f} s "
                 "(synaptic-current-weighted population signal, superficial to deep)",
                 fontsize=10, loc="left")

    # ---- (b) single-epoch zoom: one oscillation event at fine resolution ----
    ax = fig.add_subplot(gs[1])
    mid = (window[0] + window[1]) / 2.0
    ew = (bins >= mid - epoch_ms / 2) & (bins < mid + epoch_ms / 2)
    ax.plot(bins[ew], composite[ew], color="0.35", lw=1.0, label="composite cortical LFP")
    ax.plot(bins[ew], thal[ew], color=BLUE, lw=1.0, alpha=0.8, label="thalamic (TC+RE) LFP")
    ax.set_title(f"(b) Single-epoch zoom ({epoch_ms/1000:.1f} s window) -- "
                 "compare waveform shape to a literature spindle epoch",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (ms)"); ax.set_ylabel("a.u.")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    # ---- (c) composite LFP: raw + spindle band + spectrogram-style summary ----
    ax = fig.add_subplot(gs[2])
    spin = _bandpass(composite, fs, 8.0, 15.0)
    env = np.abs(hilbert(spin))
    ax.plot(tsec, composite[w], color="0.5", lw=0.6, label="composite LFP (raw)")
    ax.plot(tsec, spin[w] - 8, color=GREEN, lw=0.9, label="8-15 Hz spindle band")
    ax.plot(tsec, env[w] - 8, color="k", lw=0.8, alpha=0.6, label="envelope")
    ax.plot(tsec, -env[w] - 8, color="k", lw=0.8, alpha=0.6)
    ax.set_title("(c) Composite cortical LFP: raw vs. spindle-band-filtered + envelope",
                 fontsize=10, loc="left")
    ax.set_xlabel("time (s)"); ax.set_ylabel("a.u.")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    fig.suptitle("LFP proxy -- for direct comparison against published spindle/SO figures",
                 fontsize=13, y=0.995)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved LFP figure to {out_png}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("npz", type=str, help="path to a ctx_thalamus_mpi.py --out .npz")
    ap.add_argument("--outdir", type=str, default="out")
    ap.add_argument("--tag", type=str, default=None,
                     help="output filename prefix (default: derived from the npz filename)")
    ap.add_argument("--window-start", type=float, default=20000.0,
                     help="zoomed-window start (ms), default 20000")
    ap.add_argument("--window-len", type=float, default=10000.0,
                     help="zoomed-window length (ms), default 10000")
    args = ap.parse_args(argv)

    npz = np.load(args.npz, allow_pickle=True)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    tag = args.tag or Path(args.npz).stem
    window = (args.window_start, args.window_start + args.window_len)

    make_meanfield_figure(npz, outdir / f"{tag}_meanfield.png", window=window)
    s = make_spindle_figure(npz, outdir / f"{tag}_spindles.png", window=window)
    make_lfp_figure(npz, outdir / f"{tag}_lfp.png", window=window)

    print(f"RE burst shape: frac_burst={s['frac_burst']:.3f} "
          f"mean_burst_size={s['mean_burst_size']:.2f} n_events={s['n_events']} "
          f"event_Hz={s['event_hz']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
