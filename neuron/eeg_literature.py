"""
Sleep-EEG graphs of the model's EEG proxy (eeg_proxy.py), made the way they are
made in Moelle et al. 2002 (J Neurosci 22:10941) and in sleep-scoring displays
(30 s epochs, 1 s grid), so that they can be compared directly with the
literature. All amplitudes are z-units, so only shapes, timing and unit-free
ratios are comparable.

Moelle et al. procedure used here:
  - slow-oscillation signal: 0.16-4 Hz; spindle signal: 12-15 Hz; beta 15-25 Hz,
    alpha 8-12 Hz, theta 4-8 Hz; rms of each band in a 0.1 s window;
  - half-waves: two successive zero crossings of the SO signal 0.125-1 s apart,
    the peak between them beyond a threshold. Moelle used 80 uV and raised it per
    subject until more than 50 waves were found; the proxy has no microvolts, so
    the threshold here is the median peak amplitude of each polarity per run
    (the larger half of the half-waves);
  - wave-triggered averages in +-1 s around the peak (Fig 4);
  - spindles: 12-15 Hz rms, edges at its mean, peak above mean + 1.5 SD, 0.4-1.3 s, centre-triggered averages of the SO signal and rms (Fig 5);
  - cross-correlation of the SO signal and the spindle rms, +-2 s (Fig 6).
The head filter here uses a 0.05 Hz high-pass (the scalp study used DC recording).

NOISE (--noise R, --alpha A): an optional background with a 1/f^A power spectrum (A = 1 by
default) is added to the proxy to stand in for sources the model does not have (other
cortical areas, other generators, instrument noise). R is its rms relative to the proxy
rms in 0.05-35 Hz. R is a FREE PARAMETER: it is not measured and not fitted; the report
scans 0, 0.5, 1 and 2. With R = 2 the unit-free ratios come close to the human values,
but that is a consequence of the scan, not a validation. Noise cannot create the negative
half-wave of a scalp slow wave; only a different source (membrane or synaptic current)
could. One channel only, so no channel-common noise or topography is modelled.

    python3 neuron/eeg_literature.py --run run0.npz ... --out-prefix figures/b120_lit \
        --out-txt res/b120_lit_table.txt
"""

import argparse

import numpy as np
from scipy.signal import welch

import eeg_proxy as P

FS = P.FS
SKIP = int(P.SKIP_MS)
BANDS = {"spindle 12-15 Hz": (12, 15), "beta 15-25 Hz": (15, 25), "alpha 8-12 Hz": (8, 12),
         "theta 4-8 Hz": (4, 8)}


def rms(s, win_ms=100):
    k = int(win_ms * FS / 1000)
    return np.sqrt(np.convolve(s ** 2, np.ones(k) / k, mode="same"))


def add_noise(x, ratio, alpha=1.0, seed=0):
    """Background of non-modelled sources: noise with a 1/f^alpha power spectrum, rms =
    ratio x rms of the signal (both in the 0.05-35 Hz band). ratio is a free parameter."""
    if ratio <= 0:
        return x
    rng = np.random.default_rng(seed)
    n = len(x)
    f = np.fft.rfftfreq(n, 1 / FS)
    spec = np.fft.rfft(rng.standard_normal(n))
    g = np.zeros_like(f)
    g[1:] = f[1:] ** (-alpha / 2.0)
    nz = np.fft.irfft(spec * g, n)
    nz = P.head_filter(nz, band=(0.05, 35.0))
    return x + nz * ratio * x[SKIP:].std() / nz[SKIP:].std()


def prepare(path, use=None, noise=0.0, alpha=1.0, seed=0):
    z = np.load(path, allow_pickle=True)
    x = P.head_filter(P.source(z, use=use), band=(0.05, 35.0))
    x = add_noise(x, noise, alpha, seed)
    so = P.bandpass(x, 0.16, 4.0, 2)
    r = {name: rms(P.bandpass(x, a, b)) for name, (a, b) in BANDS.items()}
    return dict(x=x, so=so, rms=r, T=float(z["tstop"]))


def half_waves(so, sign, lo=0.125, hi=1.0):
    """Peak indices and amplitudes of positive (sign=+1) or negative (-1) half-waves."""
    y = sign * so
    zc = np.where((y[:-1] <= 0) & (y[1:] > 0))[0]       # starts of the half-waves
    ze = np.where((y[:-1] > 0) & (y[1:] <= 0))[0]
    out = []
    for a in zc:
        e = ze[ze > a]
        if len(e) == 0:
            break
        b = e[0]
        d = (b - a) / FS
        if lo <= d <= hi:
            pk = a + int(np.argmax(y[a:b]))
            out.append((pk, y[pk], a, b))
    return np.array(out).reshape(-1, 4)


def select(hw):
    """Larger half of the half-waves (>= the median peak), after the first SKIP ms."""
    hw = hw[hw[:, 0] >= SKIP]
    return hw[hw[:, 1] >= np.median(hw[:, 1])]


def triggered(sig, idx, w=1000):
    return np.array([sig[int(i) - w:int(i) + w + 1] for i in idx
                     if i - w >= 0 and i + w < len(sig)])


def spindles_moelle(r12, thr_sd=1.5, lo=0.4, hi=1.3):
    """Spindles on the 12-15 Hz rms: edges where the rms crosses its mean, peak above
    mean + thr_sd SD, duration 0.4-1.3 s (the scalp study used one threshold; with one
    threshold the short model bursts are not found, see the report)."""
    ok = np.arange(len(r12)) >= SKIP
    mu, sd = r12[ok].mean(), r12[ok].std()
    a = (r12 > mu) & ok
    d = np.diff(a.astype(int))
    on, off = np.where(d == 1)[0] + 1, np.where(d == -1)[0] + 1
    if len(off) and len(on) and off[0] < on[0]:
        off = off[1:]
    n = min(len(on), len(off))
    c = []
    for u, v in zip(on[:n], off[:n]):
        if lo <= (v - u) / FS <= hi and r12[u:v].max() >= mu + thr_sd * sd:
            c.append((u + v) // 2)
    return np.array(c, int)


def xcorr(so, rs, max_lag_s=2.0, step=50):
    """corr(SO(t), rms(t + lag)), series on a 0.05 s grid, over the run after SKIP."""
    a = so[SKIP::step]
    b = rs[SKIP::step]
    a = (a - a.mean()) / a.std()
    b = (b - b.mean()) / b.std()
    n = int(max_lag_s * FS / step)
    lags = np.arange(-n, n + 1)
    c = []
    for l in lags:
        if l >= 0:
            c.append(np.mean(a[:len(a) - l] * b[l:]))
        else:
            c.append(np.mean(a[-l:] * b[:len(b) + l]))
    return lags * step / FS, np.array(c)


def analyse(D):
    """Per-run and pooled results used by the table and the figures."""
    res = dict(neg=[], pos=[], spin_c=[], xc=[], so_psd=[], n_neg=[], n_pos=[], n_sp=[], runs=D)
    for d in D:
        for sign, key in ((-1, "neg"), (1, "pos")):
            hw = select(half_waves(d["so"], sign))
            idx = hw[:, 0].astype(int)
            res["n_" + key].append(len(idx))
            res[key].append({"so": triggered(d["so"], idx), **{n: triggered(r, idx)
                                                                for n, r in d["rms"].items()}})
        c = spindles_moelle(d["rms"]["spindle 12-15 Hz"])
        res["n_sp"].append(len(c))
        res["spin_c"].append({"so": triggered(d["so"], c),
                              "rms": triggered(d["rms"]["spindle 12-15 Hz"], c)})
        res["xc"].append(xcorr(d["so"], d["rms"]["spindle 12-15 Hz"])[1])
        f, p = welch(d["so"][SKIP:], fs=FS, nperseg=int(30 * FS), noverlap=int(15 * FS))
        res["so_psd"].append((f, p))
    res["lags"] = xcorr(D[0]["so"], D[0]["rms"]["spindle 12-15 Hz"])[0]
    return res


def _pool(lst, key):
    """Per-run means of the triggered arrays (runs are the units)."""
    return np.array([x[key].mean(0) for x in lst if len(x[key])])


def table(R):
    tt = np.arange(-1000, 1001) / 1000.0
    base = tt < -0.8
    out = ["Moelle-style analysis of the EEG proxy, %d runs of 120 s (mean +- SEM over runs)"
           % len(R["runs"]), ""]
    out.append("half-waves per run (larger half): negative %.0f +- %.0f, positive %.0f +- %.0f; "
               "spindles per run (12-15 Hz rms, 0.4-1.3 s): %.0f +- %.0f"
               % (np.mean(R["n_neg"]), np.std(R["n_neg"], ddof=1), np.mean(R["n_pos"]),
                  np.std(R["n_pos"], ddof=1), np.mean(R["n_sp"]), np.std(R["n_sp"], ddof=1)))

    def ms(v, d=2):
        v = np.asarray(v, float)
        return "%.*f +- %.*f" % (d, v.mean(), d, v.std(ddof=1) / np.sqrt(len(v)))
    # mean half-wave length: width of the average between the zero crossings around the peak
    rows = {}
    for band in BANDS:
        pn = _pool(R["neg"], band)
        pp = _pool(R["pos"], band)
        rows[band] = (pn, pp)
    sn, sp = rows["spindle 12-15 Hz"]
    i0 = 1000
    out.append("")
    out.append("Spindle rms (z-units) at the peak of the half-wave: negative %s, positive %s, "
               "ratio positive/negative %s" % (ms(sn[:, i0], 3), ms(sp[:, i0], 3),
                                               ms(sp[:, i0] / sn[:, i0])))
    out.append("  Moelle et al.: Cz 3.32 vs 2.30 uV (ratio 1.44), Fz 3.13 vs 2.26 uV (ratio 1.39)")
    b_n, b_p = sn[:, base].mean(1), sp[:, base].mean(1)
    out.append("Change of spindle rms against the first 0.2 s: negative %s %%, positive %s %%"
               % (ms(100 * (sn[:, i0] / b_n - 1), 0), ms(100 * (sp[:, i0] / b_p - 1), 0)))
    out.append("  Moelle et al.: -22 % and +23.4 % at Cz")
    w = (tt >= 0.2) & (tt <= 0.8)
    rb_t = tt[w][np.argmax(sn[:, w].mean(0))]
    out.append("Rebound after the negative half-wave: spindle rms peaks at +%.2f s, value %s against "
               "%s at the positive peak" % (rb_t, ms(sn[:, w].max(1), 3), ms(sp[:, i0], 3)))
    out.append("  Moelle et al.: rebound 400-500 ms after the trough, higher than at the positive peak")
    for band in ("beta 15-25 Hz", "alpha 8-12 Hz", "theta 4-8 Hz"):
        bn, bp = rows[band]
        out.append("%s rms at the peak: negative %s, positive %s (ratio %s)"
                   % (band, ms(bn[:, i0], 3), ms(bp[:, i0], 3), ms(bp[:, i0] / bn[:, i0])))
    # spindle-triggered
    sc = _pool(R["spin_c"], "so")
    sr = _pool(R["spin_c"], "rms")
    bl = tt < -0.8
    sc0 = sc - sc[:, bl].mean(1, keepdims=True)
    mneg = tt[(tt >= -0.8) & (tt <= 0)][np.argmin(sc0[:, (tt >= -0.8) & (tt <= 0)].mean(0))]
    mpos = tt[(tt >= -0.5) & (tt <= 0.5)][np.argmax(sc0[:, (tt >= -0.5) & (tt <= 0.5)].mean(0))]
    out.append("")
    out.append("Spindle-triggered SO signal: negative peak at %.2f s before the spindle centre, "
               "positive peak at %+.2f s" % (-mneg, mpos))
    out.append("  Moelle et al.: negative peak shortly before the rise of spindle activity, positive "
               "peak with the highest spindle power, then negative again")
    xc = np.array(R["xc"])
    lags = R["lags"]
    m = xc.mean(0)
    se = xc.std(0, ddof=1) / np.sqrt(len(xc))
    sel = (lags >= -1.0) & (lags <= 1.0)
    lneg = lags[sel][np.argmin(m[sel])]
    lpos = lags[(lags >= -0.5) & (lags <= 0.5)][np.argmax(m[(lags >= -0.5) & (lags <= 0.5)])]
    out.append("Cross-correlation SO vs spindle rms: strongest positive peak at %+.2f s (r = %.2f), "
               "strongest negative peak at %+.2f s (r = %.2f, SEM %.2f)"
               % (lpos, m[(lags >= -0.5) & (lags <= 0.5)].max(), lneg, m[sel].min(),
                  se[sel][np.argmin(m[sel])]))
    out.append("  Moelle et al.: positive peak at -0.1 s, strongest (negative) peak at +0.4 s")
    f, _ = R["so_psd"][0]
    P_ = np.array([p for _, p in R["so_psd"]])
    s = (f >= 0.16) & (f <= 4)
    pk = [float(f[s][np.argmax(p[s])]) for p in P_]
    out.append("")
    out.append("Spectral peak of the 0.16-4 Hz signal (30 s Welch): %s Hz" % ms(pk, 2))
    out.append("  Moelle et al.: single broad peak at 0.7-0.8 Hz, small one at about 0.25 Hz")
    return out


def figures(R, pre):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    D = R["runs"]
    tt = np.arange(-1000, 1001) / 1000.0

    # ---- epoch display (cf. Moelle Fig 1/3 and the N3 slides): 30 s, 1 s grid ----
    d = D[0]
    t0 = 30
    s = slice(t0 * 1000, (t0 + 30) * 1000)
    t = np.arange(len(d["x"]))[s] / FS
    pos = half_waves(d["so"], 1)
    neg = half_waves(d["so"], -1)
    posS, negS = select(pos), select(neg)
    fig, ax = plt.subplots(5, 1, figsize=(15, 9), sharex=True,
                           gridspec_kw={"height_ratios": [1.4, 1.2, 1, 1, 1]})
    for a in ax:
        for g in range(t0, t0 + 31):
            a.axvline(g, color="0.85", lw=0.5, zorder=0)
    ax[0].plot(t, d["x"][s], "k", lw=0.6)
    ax[0].set_ylabel("EEG proxy\n0.05-35 Hz")
    ax[1].plot(t, d["so"][s], "k", lw=0.9)
    for hw, c in ((posS, "C3"), (negS, "C0")):
        for pk, _, a_, b_ in hw:
            if t0 * 1000 <= a_ and b_ < (t0 + 30) * 1000:
                ax[1].plot(np.arange(int(a_), int(b_) + 1) / FS, d["so"][int(a_):int(b_) + 1],
                           c, lw=2.2)
    ax[1].set_ylabel("slow oscillation\n0.16-4 Hz")
    sp_band = P.bandpass(d["x"], 12, 15)
    ax[2].plot(t, sp_band[s], "C2", lw=0.6)
    ax[2].set_ylabel("spindle\n12-15 Hz")
    ax[3].plot(t, d["rms"]["spindle 12-15 Hz"][s], "C2")
    cs = spindles_moelle(d["rms"]["spindle 12-15 Hz"])
    for c in cs:
        if t0 * 1000 <= c < (t0 + 30) * 1000:
            ax[3].plot(c / FS, d["rms"]["spindle 12-15 Hz"][c], "v", color="C1", ms=6)
    ax[3].set_ylabel("spindle rms\n(0.1 s)")
    ax[4].plot(t, d["rms"]["beta 15-25 Hz"][s], "C4")
    ax[4].set_ylabel("beta rms\n15-25 Hz")
    ax[4].set_xlabel("time (s), grid = 1 s")
    ax[0].set_title("30 s epoch (cf. Moelle Fig 1/3, sleep-scoring displays): red/blue = larger positive/"
                    "negative half-waves, orange = spindle centres")
    fig.tight_layout()
    fig.savefig(pre + "_epoch.png", dpi=120)
    plt.close(fig)

    # ---- Fig 2: spectrum of the SO signal ----
    f = R["so_psd"][0][0]
    Pm = np.array([p for _, p in R["so_psd"]])
    s = (f >= 0.05) & (f <= 4.5)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(f[s], Pm.mean(0)[s], "k")
    ax.fill_between(f[s], (Pm.mean(0) - Pm.std(0, ddof=1))[s], (Pm.mean(0) + Pm.std(0, ddof=1))[s],
                    color="k", alpha=0.2)
    ax.axvspan(0.7, 0.8, color="C1", alpha=0.3, label="Moelle et al. peak 0.7-0.8 Hz")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("power density (z^2/Hz)")
    ax.set_title("Spectrum of the 0.16-4 Hz signal (cf. Moelle Fig 2), %d runs" % len(D))
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(pre + "_spectrum.png", dpi=130)
    plt.close(fig)

    # ---- Fig 4: wave-triggered averages ----
    keys = ["so", "spindle 12-15 Hz", "beta 15-25 Hz", "alpha 8-12 Hz", "theta 4-8 Hz"]
    ylab = ["slow oscillation (z)", "spindle rms 12-15 Hz (z)", "beta rms 15-25 Hz (z)",
            "alpha rms 8-12 Hz (z)", "theta rms 4-8 Hz (z)"]
    fig, ax = plt.subplots(5, 2, figsize=(10, 12), sharex=True)
    for col, (lst, ttl) in enumerate(((R["neg"], "negative half-waves"),
                                      (R["pos"], "positive half-waves"))):
        for row, k in enumerate(keys):
            X = _pool(lst, k)
            m, e = X.mean(0), X.std(0, ddof=1) / np.sqrt(len(X))
            ax[row, col].plot(tt, m, "k")
            ax[row, col].fill_between(tt, m - e, m + e, color="k", alpha=0.25)
            ax[row, col].axvline(0, color="C3", ls=":", lw=0.8)
            if col == 0:
                ax[row, col].set_ylabel(ylab[row])
        ax[0, col].set_title("%s (%d per run)" % (ttl, np.mean(R["n_neg" if col == 0 else "n_pos"])))
        ax[4, col].set_xlabel("time from the half-wave peak (s)")
    fig.suptitle("Half-wave-triggered averages, mean +- SEM over %d runs (cf. Moelle Fig 4)" % len(D))
    fig.tight_layout()
    fig.savefig(pre + "_halfwave_avg.png", dpi=120)
    plt.close(fig)

    # ---- Fig 5: spindle-triggered ----
    sc = _pool(R["spin_c"], "so")
    sr = _pool(R["spin_c"], "rms")
    sc = sc - sc[:, tt < -0.8].mean(1, keepdims=True)
    fig, ax = plt.subplots(2, 1, figsize=(6.5, 6), sharex=True)
    for a, X, yl in ((ax[0], sr, "spindle rms (z)"), (ax[1], sc, "SO signal (z), baseline-corrected")):
        m, e = X.mean(0), X.std(0, ddof=1) / np.sqrt(len(X))
        a.plot(tt, m, "k")
        a.fill_between(tt, m - e, m + e, color="k", alpha=0.25)
        a.axvline(0, color="C3", ls=":", lw=0.8)
        a.set_ylabel(yl)
    ax[1].set_xlabel("time from the spindle centre (s)")
    ax[0].set_title("Spindle-triggered averages, %d runs (cf. Moelle Fig 5)" % len(D))
    fig.tight_layout()
    fig.savefig(pre + "_spindle_avg.png", dpi=130)
    plt.close(fig)

    # ---- Fig 6: cross-correlation ----
    xc = np.array(R["xc"])
    lags = R["lags"]
    fig, ax = plt.subplots(figsize=(6.5, 4))
    m, e = xc.mean(0), xc.std(0, ddof=1) / np.sqrt(len(xc))
    ax.plot(lags, m, "k")
    ax.fill_between(lags, m - e, m + e, color="k", alpha=0.25)
    ax.axvline(0, color="C3", ls=":", lw=0.8)
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlabel("lag (s); positive = spindle rms lags the SO signal")
    ax.set_ylabel("cross-correlation")
    ax.set_title("SO signal vs spindle rms, %d runs (cf. Moelle Fig 6)" % len(D))
    fig.tight_layout()
    fig.savefig(pre + "_xcorr.png", dpi=130)
    plt.close(fig)


def noise_figure(path, use, pre, levels=(0.0, 1.0, 2.0), t0=30, alpha=1.0):
    """Same 30 s of run 0 at several noise levels, and the spectra, to compare with a
    real N3 recording."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(15, 8))
    gs = fig.add_gridspec(len(levels), 2, width_ratios=[2.6, 1])
    spec = []
    for i, lv in enumerate(levels):
        d = prepare(path, use, lv, alpha, seed=0)
        a = fig.add_subplot(gs[i, 0])
        s = slice(t0 * 1000, (t0 + 30) * 1000)
        a.plot(np.arange(len(d["x"]))[s] / FS, d["x"][s], "k", lw=0.6)
        for g in range(t0, t0 + 31):
            a.axvline(g, color="0.85", lw=0.5, zorder=0)
        a.set_ylabel("noise rms = %g x signal" % lv)
        if i == 0:
            a.set_title("EEG proxy (superficial E cells) plus 1/f background, same 30 s of run 0")
        f, pw = welch(d["x"][SKIP:], fs=FS, nperseg=int(8 * FS), noverlap=int(4 * FS))
        spec.append((lv, f, pw))
    a.set_xlabel("time (s), grid = 1 s")
    b = fig.add_subplot(gs[:, 1])
    for lv, f, pw in spec:
        m = (f >= 0.3) & (f <= 40)
        b.loglog(f[m], pw[m], label="noise %g" % lv)
    b.set_xlabel("frequency (Hz)")
    b.set_ylabel("power (z^2/Hz)")
    b.legend(frameon=False)
    b.set_title("spectra")
    fig.tight_layout()
    fig.savefig(pre + "_noise.png", dpi=120)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--out-txt")
    ap.add_argument("--noise", type=float, default=0.0,
                    help="rms of the 1/f background relative to the proxy (default 0)")
    ap.add_argument("--noise-figure", action="store_true", help="also write the noise comparison")
    ap.add_argument("--alpha", type=float, default=1.0, help="spectral exponent of the noise")
    ap.add_argument("--layers", help="comma list of E groups, e.g. l23e,l4e (default all four)")
    a = ap.parse_args(argv)
    use = a.layers.split(",") if a.layers else None
    R = analyse([prepare(p, use, a.noise, a.alpha, seed=i) for i, p in enumerate(a.run)])
    txt = table(R)
    print("\n".join(txt))
    if a.out_txt:
        open(a.out_txt, "w").write("\n".join(txt) + "\n")
    figures(R, a.out_prefix)
    if a.noise_figure:
        noise_figure(a.run[0], use, a.out_prefix, alpha=a.alpha)


if __name__ == "__main__":
    main()
