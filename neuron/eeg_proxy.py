"""
Scalp-EEG proxy of the undisturbed model and the standard sleep-EEG graphs
computed from it, for comparison with the literature.

The model has no geometry, so a head model is not possible. The proxy is a
deliberately simple chain:
  1. source: spikes of the cortical excitatory cells (pyramidal cells dominate
     scalp EEG) of L2/3, L4, L5, L6, each convolved with a PSP kernel
     (2 ms rise, 20 ms decay by default; the 100 ms kernel of lfp_baseline.py
     removes about 8/10 of the 12 Hz amplitude), z-scored, L5 and L6 sign-flipped
     as in lfp_baseline.py; the EEG source is their mean;
  2. skull and scalp: a first-order low-pass (default corner 40 Hz) that stands
     for the spatial smoothing and the attenuation of fast activity, then the
     usual recording band 0.3-35 Hz (4th order, zero phase). Averaging over the
     whole population already plays the role of spatial smoothing over a large
     cortical patch;
  3. units are z-units (SD of the unfiltered source), not microvolts.
Not modelled: the head geometry, source orientation, other generators (hippocampus,
non-modelled cortex), the measurement noise.

Analyses, per run and pooled over runs (mean +- SD over runs):
  - Welch spectrum of the EEG proxy, band powers;
  - slow oscillations (SO): 0.16-1.25 Hz band, a SO is a negative-to-negative
    zero-crossing cycle of 0.8-2.0 s with the trough in the lowest quartile
    (Staresina et al. 2015 criteria); density, trough-to-peak amplitude;
  - SO-triggered average of the EEG and of the sigma (11-16 Hz) envelope;
  - spindles: 11-16 Hz band, Hilbert envelope smoothed over 200 ms, peak above
    mean + 1.5 SD, edges where the envelope crosses its mean, 0.5-3 s; density, duration, frequency, SO phase of the
    spindle maximum, fraction of spindles that follow a SO trough within 1.5 s.

    python3 neuron/eeg_proxy.py --run run0.npz --run run1.npz ... \
        --out-prefix figures/b120_eeg --out-txt res/b120_eeg_table.txt
"""

import argparse

import numpy as np
from scipy.signal import welch, butter, sosfiltfilt, sosfilt, hilbert, find_peaks

import ctx_analyze as A
import lfp_baseline as L

E_LAYERS = [("L2/3 E", "l23e", False), ("L4 E", "l4e", False), ("L5 E", "l5e", True),
            ("L6 E", "l6e", True)]
FS = 1000.0
SKIP_MS = 2000.0


def source(z, tau_decay=20.0, use=None):
    """Mean of the z-scored, sign-corrected PSP traces of the four E populations."""
    t, g, R, T = z["times"], z["gids"], z["ranges"].item(), float(z["tstop"])
    xs = []
    for _, k, flip in E_LAYERS:
        if use is not None and k not in use:
            continue
        x, _ = A._synaptic_lfp(L._pool(t, g, R, [k])[0], T, tau_decay=tau_decay,
                               kernel_ms=10 * tau_decay)
        x = (x - x.mean()) / (x.std() + 1e-9)
        xs.append(-x if flip else x)
    return np.mean(xs, axis=0)


def head_filter(x, skull_hz=40.0, band=(0.3, 35.0)):
    """Skull/scalp low-pass followed by the recording band."""
    sos = butter(1, skull_hz, "low", fs=FS, output="sos")
    y = sosfilt(sos, x)
    sos = butter(4, band, "band", fs=FS, output="sos")
    return sosfiltfilt(sos, y)


def bandpass(x, lo, hi, order=3):
    return sosfiltfilt(butter(order, [lo, hi], "band", fs=FS, output="sos"), x)


def detect_so(x):
    """Staresina-like SO detection on the 0.16-1.25 Hz band. Returns trough
    indices, trough-to-next-peak amplitudes, cycle durations (s)."""
    s = bandpass(x, 0.16, 1.25, 2)
    neg = np.where((s[:-1] > 0) & (s[1:] <= 0))[0]      # positive-to-negative crossings
    cyc = []
    for a, b in zip(neg[:-1], neg[1:]):
        d = (b - a) / FS
        if 0.8 <= d <= 2.0:
            seg = s[a:b]
            tr = a + int(np.argmin(seg))
            pk = a + int(np.argmax(seg))
            cyc.append((tr, seg.max() - seg.min(), d, seg.min(), pk))
    if not cyc:
        return np.array([], int), np.array([]), np.array([]), s, np.array([], int), 0
    cyc = np.array(cyc)
    n_all = len(cyc)
    keep = cyc[:, 3] <= np.percentile(cyc[:, 3], 25)    # trough in the lowest quartile
    return cyc[keep, 0].astype(int), cyc[keep, 1], cyc[keep, 2], s, cyc[keep, 4].astype(int), n_all


def detect_spindles(x, so_s):
    """Amplitude-based spindle detection in the 11-16 Hz band."""
    s = bandpass(x, 11.0, 16.0)
    env = np.abs(hilbert(s))
    k = int(0.2 * FS)
    env_s = np.convolve(env, np.ones(k) / k, mode="same")
    ok = np.arange(len(x)) >= SKIP_MS
    mu, sd = env_s[ok].mean(), env_s[ok].std()
    thr = mu + 1.5 * sd                                 # peak criterion
    above = (env_s > mu) & ok                           # edges: envelope above its mean
    d = np.diff(above.astype(int))
    on, off = np.where(d == 1)[0] + 1, np.where(d == -1)[0] + 1
    if len(off) and len(on) and off[0] < on[0]:
        off = off[1:]
    n = min(len(on), len(off))
    mon, moff = list(on[:n]), list(off[:n])
    ph = np.angle(hilbert(so_s))
    res = []
    for a, b in zip(mon, moff):
        dur = (b - a) / FS
        if not 0.5 <= dur <= 3.0:
            continue
        pk = a + int(np.argmax(env_s[a:b]))
        if env_s[pk] < thr:
            continue
        seg = x[a:b] - x[a:b].mean()
        fr = np.fft.rfftfreq(8192, 1 / FS)
        sp_ = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), 8192))
        m_ = (fr >= 9) & (fr <= 20)
        res.append((a, b, dur, fr[m_][np.argmax(sp_[m_])], ph[pk], env_s[a:b].max(), pk))
    return np.array(res).reshape(-1, 7), s, env_s, thr


def run_stats(path):
    z = np.load(path, allow_pickle=True)
    T = float(z["tstop"])
    x = head_filter(source(z))
    so_idx, so_amp, so_dur, so_s, so_pk, n_all = detect_so(x)
    keep = so_idx >= SKIP_MS
    so_idx, so_amp, so_dur, so_pk = so_idx[keep], so_amp[keep], so_dur[keep], so_pk[keep]
    sp, sig, env, thr = detect_spindles(x, so_s)
    span_min = (T - SKIP_MS) / 60000.0
    # spindle follows SO trough within 1.5 s
    follow = [bool(np.any((p[0] - so_idx >= 0) & (p[0] - so_idx <= 1.5 * FS))) for p in sp]
    f, p = welch(x[int(SKIP_MS):], fs=FS, nperseg=int(8 * FS), noverlap=int(4 * FS))
    return dict(T=T, x=x, so_s=so_s, sig=sig, env=env, thr=thr, so_idx=so_idx, so_pk=so_pk, so_all_per_min=n_all / ((T - 0) / 60000.0), so_amp=so_amp,
                so_dur=so_dur, sp=sp, follow=np.array(follow, bool), f=f, p=p,
                so_per_min=len(so_idx) / span_min, sp_per_min=len(sp) / span_min)


def collect(paths):
    return [run_stats(p) for p in paths]


def circ(ph):
    if len(ph) == 0:
        return np.nan, np.nan
    v = np.exp(1j * ph).mean()
    return float(np.degrees(np.angle(v))), float(abs(v))


def tables(D):
    def ms(v, d=2):
        v = np.asarray(v, float)
        return f"{v.mean():.{d}f} +- {v.std(ddof=1):.{d}f}"
    out = ["EEG proxy of the undisturbed model, %d runs (mean +- SD over runs)" % len(D), ""]
    out.append("SO (0.16-1.25 Hz, Staresina-like criteria): %s /min, trough-to-peak %s z-units, "
               "cycle %s s" % (ms([d["so_per_min"] for d in D]),
                               ms([d["so_amp"].mean() for d in D]),
                               ms([d["so_dur"].mean() for d in D])))
    out.append("All 0.8-2 s negative-to-negative SO cycles (before the quartile criterion): %s /min"
               % ms([d["so_all_per_min"] for d in D], 1))
    out.append("Spindles (11-16 Hz, peak > mean+1.5 SD, edges at the mean, 0.5-3 s): %s /min, duration %s s, "
               "frequency (FFT peak in 9-20 Hz) %s Hz" % (ms([d["sp_per_min"] for d in D], 1),
                                    ms([d["sp"][:, 2].mean() for d in D]),
                                    ms([d["sp"][:, 3].mean() for d in D], 1)))
    ph = np.concatenate([d["sp"][:, 4] for d in D])
    m, r = circ(ph)
    n = len(ph)
    pr = float(np.exp(-n * r * r)) if n else np.nan
    out.append("Spindle maximum at SO phase %.0f deg (0 = positive SO peak), R = %.2f, n = %d, "
               "Rayleigh p ~ %.2g" % (m, r, n, pr))
    fol = np.concatenate([d["follow"] for d in D])
    out.append("Spindles starting 0-1.5 s after an SO trough: %.0f %%" % (100 * fol.mean()))
    out.append("SOs followed by a spindle within 1.5 s: %.0f %%" % (100 * np.mean(
        [np.any((d["sp"][:, 0] - t >= 0) & (d["sp"][:, 0] - t <= 1.5 * FS)) for d in D
         for t in d["so_idx"]])))
    out.append("")
    out.append("Relative power of the EEG proxy, 0.3-35 Hz (mean +- SD over runs):")
    bands = [("SO 0.3-1.5", 0.3, 1.5), ("delta 1.5-4", 1.5, 4), ("theta 4-8", 4, 8),
             ("alpha 8-11", 8, 11), ("sigma 11-16", 11, 16), ("beta 16-30", 16, 30)]
    for name, a, b in bands:
        fr = []
        for d in D:
            f, p = d["f"], d["p"]
            tot = np.trapezoid(p[(f >= 0.3) & (f <= 35)], f[(f >= 0.3) & (f <= 35)])
            s = (f >= a) & (f < b)
            fr.append(np.trapezoid(p[s], f[s]) / tot)
        out.append("  %-12s %s" % (name, ms(fr, 3)))
    return out


def figures(D, pre):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- 1: raw trace, SO band, sigma band, detected spindles ----
    d = D[0]
    t0, t1 = 20, 50
    t = np.arange(len(d["x"])) / FS
    s = slice(int(t0 * FS), int(t1 * FS))
    fig, ax = plt.subplots(3, 1, figsize=(13, 7), sharex=True)
    ax[0].plot(t[s], d["x"][s], "k", lw=0.6)
    ax[0].set_ylabel("EEG proxy (z)")
    ax[1].plot(t[s], d["so_s"][s], "C0", lw=1)
    for i in d["so_idx"]:
        if t0 * FS <= i < t1 * FS:
            ax[1].plot(i / FS, d["so_s"][i], "v", color="C3", ms=6)
    ax[1].set_ylabel("SO band 0.16-1.25 Hz")
    ax[2].plot(t[s], d["sig"][s], "C2", lw=0.6)
    ax[2].plot(t[s], d["env"][s], "k", lw=0.8)
    ax[2].axhline(d["thr"], color="C3", ls="--", lw=0.8)
    for a, b, *_ in d["sp"]:
        if t0 * FS <= a < t1 * FS:
            ax[2].axvspan(a / FS, b / FS, color="C1", alpha=0.3)
    ax[2].set_ylabel("sigma 11-16 Hz")
    ax[2].set_xlabel("time (s)")
    ax[0].set_title("EEG proxy of run 0: raw, SO band with detected SO troughs, sigma band with "
                    "detected spindles (orange)")
    fig.tight_layout()
    fig.savefig(pre + "_traces.png", dpi=130)
    plt.close(fig)

    # ---- 2: spectrum and spectrogram ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1, 1.3]})
    f = D[0]["f"]
    P = np.array([d["p"] for d in D])
    sel = (f >= 0.3) & (f <= 40)
    mu, sd = P.mean(0), P.std(0, ddof=1)
    ax[0].loglog(f[sel], mu[sel], "k")
    ax[0].fill_between(f[sel], np.maximum(mu - sd, 1e-12)[sel], (mu + sd)[sel], color="k",
                       alpha=0.2)
    ax[0].set_xlabel("frequency (Hz)")
    ax[0].set_ylabel("power (z^2/Hz)")
    ax[0].set_title("Welch spectrum, mean +- SD over %d runs" % len(D))
    from scipy.signal import spectrogram
    fs_, ts_, Sx = spectrogram(D[0]["x"], fs=FS, nperseg=int(2 * FS), noverlap=int(1.75 * FS))
    m = fs_ <= 25
    ax[1].pcolormesh(ts_, fs_[m], 10 * np.log10(Sx[m] + 1e-12), shading="auto", cmap="viridis")
    ax[1].set_xlabel("time (s)")
    ax[1].set_ylabel("frequency (Hz)")
    ax[1].set_title("Spectrogram, run 0 (dB)")
    fig.tight_layout()
    fig.savefig(pre + "_spectra.png", dpi=130)
    plt.close(fig)

    # ---- 3: SO-triggered average of EEG and sigma envelope ----
    w = int(2.5 * FS)
    tr = np.arange(-w, w + 1) / FS
    fig, ax = plt.subplots(2, 2, figsize=(12, 6), sharex=True)
    for col, (key, lab) in enumerate((("so_idx", "SO trough"), ("so_pk", "SO positive peak"))):
        A_, E_ = [], []
        for d in D:
            for i in d[key]:
                if i - w >= 0 and i + w < len(d["x"]):
                    A_.append(d["x"][i - w:i + w + 1])
                    E_.append(d["env"][i - w:i + w + 1])
        A_, E_ = np.array(A_), np.array(E_)
        for a, X, c, yl in ((ax[0, col], A_, "k", "EEG proxy (z)"),
                            (ax[1, col], E_, "C2", "sigma envelope (z)")):
            m, e = X.mean(0), X.std(0, ddof=1) / np.sqrt(len(X))
            a.plot(tr, m, c)
            a.fill_between(tr, m - 1.96 * e, m + 1.96 * e, color=c, alpha=0.25)
            a.axvline(0, color="C3", lw=0.8, ls=":")
            a.set_ylabel(yl)
        ax[1, col].set_xlabel("time from %s (s)" % lab)
        ax[0, col].set_title("%d SOs from %d runs (95%% CI)" % (len(A_), len(D)))
    fig.suptitle("SO-triggered averages of the EEG proxy and the sigma envelope")
    fig.tight_layout()
    fig.savefig(pre + "_so_average.png", dpi=130)
    plt.close(fig)

    # ---- 4: spindle properties ----
    sp = np.concatenate([d["sp"] for d in D])
    fig, ax = plt.subplots(1, 4, figsize=(15, 3.6))
    ax[0].hist(sp[:, 2], bins=np.arange(0.5, 3.01, 0.1), color="C2")
    ax[0].set_xlabel("duration (s)")
    ax[1].hist(sp[:, 3], bins=np.arange(10, 18.1, 0.5), color="C2")
    ax[1].set_xlabel("frequency (Hz)")
    ph = np.degrees(sp[:, 4])
    ax[2].hist(ph, bins=np.arange(-180, 181, 30), color="C0")
    ax[2].set_xlabel("SO phase of the spindle maximum (deg, 0 = SO peak)")
    ax[3].bar(range(len(D)), [d["sp_per_min"] for d in D], color="C2")
    ax[3].set_xlabel("run")
    ax[3].set_ylabel("spindles / min")
    ax[0].set_ylabel("count")
    fig.suptitle("Detected spindles in the EEG proxy, %d spindles from %d runs" % (len(sp), len(D)))
    fig.tight_layout()
    fig.savefig(pre + "_spindles.png", dpi=130)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--out-txt")
    a = ap.parse_args(argv)
    D = collect(a.run)
    txt = tables(D)
    print("\n".join(txt))
    if a.out_txt:
        open(a.out_txt, "w").write("\n".join(txt) + "\n")
    figures(D, a.out_prefix)


if __name__ == "__main__":
    main()
