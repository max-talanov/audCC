"""
Field-potential proxies and activity statistics of the undisturbed model (no
auditory input), per population, per layer and cumulative.

LFP proxy per trace = spikes convolved with a PSP kernel
(ctx_analyze._synaptic_lfp), z-scored, deep layers and RE sign-flipped as in
ctx_analyze._composite_lfp. Traces:
  - 10 neuron groups: tc, re, l4e, l4i, l23e, l23i, l5e, l5i, l6e, l6i
  - 4 layers (E + I spikes pooled): L2/3, L4, L5, L6
  - cumulative: cortical composite (mean of the 4 layers), thalamic composite
    (mean of TC and RE) and the whole network (mean of the two).
The thalamic traces and the sign convention are proxies, not measured fields.

Per trace: Welch spectrum (8 s segments), peak frequency of the SO/delta range,
band fractions of the 0.3-40 Hz power. The sigma (9-16 Hz) peak and band come
from the population rate in 5 ms bins instead, because the 100 ms PSP kernel
attenuates 12 Hz about 8-fold. Per group: mean rate, rate in UP vs DOWN
windows, and the latency of the UP onset (half-height of the onset-locked
rate). Spindles are RE volley trains of >= 6 volleys (volley_stats.py), their
SO phase uses the zero-phase SO band of the cortical E rate (stim_phase.py).

    python3 neuron/lfp_baseline.py --run run0.npz --run run1.npz ... \
        --out-prefix figures/b120 --out-txt res/b120_table.txt
"""

import argparse

import numpy as np
from scipy.signal import welch, find_peaks, butter, sosfiltfilt, hilbert

import ctx_analyze as A
import stage2_report as S
import stim_phase as SPH
import volley_stats as V

GROUPS = ["l23e", "l23i", "l4e", "l4i", "l5e", "l5i", "l6e", "l6i", "tc", "re"]
LABEL = {"l23e": "L2/3 E", "l23i": "L2/3 I", "l4e": "L4 E", "l4i": "L4 I", "l5e": "L5 E",
         "l5i": "L5 I", "l6e": "L6 E", "l6i": "L6 I", "tc": "TC", "re": "RE"}
LAYERS = [("L2/3", ("l23e", "l23i")), ("L4", ("l4e", "l4i")), ("L5", ("l5e", "l5i")),
          ("L6", ("l6e", "l6i"))]
FLIP = {"l5e", "l5i", "l6e", "l6i", "re"}
BANDS = [("SO 0.3-1.5", 0.3, 1.5), ("delta 1.5-4", 1.5, 4.0), ("theta 4-8", 4.0, 8.0),
         ("sigma 9-16", 9.0, 16.0), ("beta 16-40", 16.0, 40.0)]
SKIP = 2000.0
UP_WIN, DOWN_WIN = (0.0, 300.0), (-400.0, -100.0)


def _pool(t, g, R, names):
    m = np.zeros_like(t, bool)
    n = 0
    for k in names:
        lo, hi = R[k]
        m |= (g >= lo) & (g < hi)
        n += hi - lo
    return t[m], n


def lfp_traces(z):
    """name -> z-scored LFP proxy at 1 kHz (flipped where needed), skip not applied."""
    t, g, R, T = z["times"], z["gids"], z["ranges"].item(), float(z["tstop"])
    out = {}

    def one(times, flip):
        x, _ = A._synaptic_lfp(times, T)
        x = (x - x.mean()) / (x.std() + 1e-9)
        return -x if flip else x

    for k in GROUPS:
        out[k] = one(*_pool(t, g, R, [k])[:1], k in FLIP)
    for name, ks in LAYERS:
        out[name] = one(_pool(t, g, R, ks)[0], ks[0] in FLIP)
    out["cortex"] = np.mean([out[n] for n, _ in LAYERS], axis=0)
    out["thalamus"] = np.mean([out["tc"], out["re"]], axis=0)
    out["network"] = np.mean([out["cortex"], out["thalamus"]], axis=0)
    return out


def psd(x, fs, nperseg_s=8.0):
    f, p = welch(x - x.mean(), fs=fs, nperseg=int(nperseg_s * fs), noverlap=int(nperseg_s * fs / 2))
    return f, p


def band_fractions(f, p, lo=0.3, hi=40.0):
    sel = (f >= lo) & (f <= hi)
    tot = np.trapezoid(p[sel], f[sel])
    out = {}
    for name, a, b in BANDS:
        s = (f >= a) & (f < b)
        out[name] = float(np.trapezoid(p[s], f[s]) / tot)
    return out


def peak_freq(f, p, lo, hi):
    s = (f >= lo) & (f <= hi)
    return float(f[s][np.argmax(p[s])])


def rate_trace(t, n, T, bin_ms):
    bins = np.arange(0, T + bin_ms, bin_ms)
    return np.histogram(t, bins)[0] / n / (bin_ms / 1000.0)


def group_stats(z, ups):
    """Mean rate, UP/DOWN rates and UP-onset profile per group."""
    t, g, R, T = z["times"], z["gids"], z["ranges"].item(), float(z["tstop"])
    prof_t = np.arange(-600, 801, 10.0)
    out = {}
    for k in GROUPS:
        tk, n = _pool(t, g, R, [k])
        tk = tk[tk >= SKIP]
        mean = len(tk) / n / ((T - SKIP) / 1000.0)
        # onset-locked rate, 10 ms bins
        prof = np.zeros(len(prof_t))
        ok = 0
        for u in ups:
            if u + prof_t[0] < SKIP or u + prof_t[-1] + 10 > T:
                continue
            h = np.histogram(tk, np.append(prof_t, prof_t[-1] + 10.0) + u)[0]
            prof += h
            ok += 1
        prof = prof / max(ok, 1) / n / 0.010
        up = prof[(prof_t >= UP_WIN[0]) & (prof_t < UP_WIN[1])].mean()
        dn = prof[(prof_t >= DOWN_WIN[0]) & (prof_t < DOWN_WIN[1])].mean()
        # half-height latency of the rise from the DOWN level to the peak in [0, 300] ms
        s = (prof_t >= -200) & (prof_t <= 300)
        base, pk = dn, prof[s].max()
        half = base + 0.5 * (pk - base)
        idx = np.nonzero(prof[s] >= half)[0]
        lat = float(prof_t[s][idx[0]]) if len(idx) and pk > base else np.nan
        out[k] = {"mean": mean, "up": up, "down": dn, "latency": lat, "profile": prof}
    return prof_t, out


def spindle_stats(z, T):
    R = z["ranges"].item()
    tr = V.volley_trains(V.re_volleys(z["times"], z["gids"], R, T, skip_ms=SKIP))
    sp = [x for x in tr if x["cycles"] >= V.SPINDLE_MIN_CYCLES]
    dur = (T - SKIP) / 1000.0
    ph, _ = SPH.phase_acausal(SPH.cortex_rate(z))
    pt = np.array([ph[int(x["t0"] // SPH.BIN_MS)] for x in sp]) if sp else np.empty(0)
    return {"n": len(sp), "per_min": len(sp) / dur * 60.0,
            "dur_s": float(np.median([x["duration"] for x in sp]) / 1000.0) if sp else np.nan,
            "cycles": float(np.median([x["cycles"] for x in sp])) if sp else np.nan,
            "freq": float(np.median([x["freq"] for x in sp])) if sp else np.nan,
            "phases": pt}


def run_stats(path):
    z = np.load(path, allow_pickle=True)
    T = float(z["tstop"])
    t, g, R = z["times"], z["gids"], z["ranges"].item()
    ups = S.up_onsets(t, g, R, T, SKIP)
    iu = np.diff(ups) / 1000.0
    L = lfp_traces(z)
    i0 = int(SKIP)
    spec = {}
    for name, x in L.items():
        f, p = psd(x[i0:], 1000.0)
        spec[name] = (f, p)
    # sigma band from the rate (5 ms bins, 200 Hz)
    rate_spec = {}
    for name in GROUPS + [n for n, _ in LAYERS]:
        ks = [name] if name in GROUPS else dict(LAYERS)[name]
        tk, n = _pool(t, g, R, ks)
        r = rate_trace(tk, n, T, 5.0)[int(SKIP / 5.0):]
        rate_spec[name] = psd(r, 200.0)
    prof_t, gs = group_stats(z, ups)
    return {"T": T, "ups": len(ups), "up_rate": len(ups) / ((T - SKIP) / 1000.0),
            "cv": float(iu.std() / iu.mean()), "spec": spec, "rate_spec": rate_spec,
            "groups": gs, "prof_t": prof_t, "spindles": spindle_stats(z, T),
            "trace": {k: L[k][20000:32000] for k in L}}


def collect(paths):
    return [run_stats(p) for p in paths]


def _ms(v, d=2):
    v = np.asarray(v, float)
    return f"{np.nanmean(v):.{d}f} ± {np.nanstd(v):.{d}f}"


def tables(D):
    lines = []
    lines.append(f"{len(D)} runs, {D[0]['T']/1000:.0f} s each, first {SKIP/1000:.0f} s skipped; mean ± SD over runs")
    lines.append("")
    lines.append("Slow oscillation and spindles")
    lines.append(f"  UP onsets per s      {_ms([d['up_rate'] for d in D])}")
    lines.append(f"  CV of UP intervals   {_ms([d['cv'] for d in D])}")
    lines.append(f"  spindles per min     {_ms([d['spindles']['per_min'] for d in D], 1)}")
    lines.append(f"  spindle duration s   {_ms([d['spindles']['dur_s'] for d in D])}")
    lines.append(f"  cycles per spindle   {_ms([d['spindles']['cycles'] for d in D], 1)}")
    lines.append(f"  spindle freq Hz      {_ms([d['spindles']['freq'] for d in D], 1)}")
    ph = np.concatenate([d["spindles"]["phases"] for d in D])
    z = np.exp(1j * ph).mean()
    lines.append(f"  spindle start SO phase: n={len(ph)}, mean {np.degrees(np.angle(z)):.0f} deg "
                 f"(0 = UP maximum), R={abs(z):.2f}, p={np.exp(-len(ph) * abs(z)**2):.2g} (Rayleigh)")
    lines.append("")
    lines.append("Per neuron group (rate Hz/cell; UP = 0-300 ms after UP onset, DOWN = -400..-100 ms; latency = half rise, ms)")
    lines.append(f"{'group':>8} | {'mean rate':>14} | {'UP rate':>14} | {'DOWN rate':>14} | {'UP/DOWN':>8} | {'latency':>14}")
    for k in GROUPS:
        m = [d["groups"][k] for d in D]
        ratio = np.mean([x["up"] for x in m]) / max(np.mean([x["down"] for x in m]), 1e-9)
        lines.append(f"{LABEL[k]:>8} | {_ms([x['mean'] for x in m]):>14} | {_ms([x['up'] for x in m]):>14} | "
                     f"{_ms([x['down'] for x in m]):>14} | {ratio:8.1f} | {_ms([x['latency'] for x in m], 0):>14}")
    lines.append("")
    lines.append("LFP proxy spectra: fraction of 0.3-40 Hz power (mean ± SD over runs)")
    names = [n for n, _ in LAYERS] + ["cortex", "thalamus", "network"] + GROUPS
    hdr = f"{'trace':>10} | " + " | ".join(f"{b[0]:>11}" for b in BANDS)
    lines.append(hdr)
    for name in names:
        fr = {b[0]: [] for b in BANDS}
        for d in D:
            f, p = d["spec"][name]
            for b, v in band_fractions(f, p).items():
                fr[b].append(v)
        lines.append(f"{LABEL.get(name, name):>10} | " +
                     " | ".join(f"{_ms(fr[b[0]]):>11}" for b in BANDS))
    lines.append("")
    lines.append("Population-rate spectra (5 ms bins): highest peak between 5 and 30 Hz (Hz) and sigma-band (9-16 Hz) fraction of 0.3-40 Hz power")
    for name in [n for n, _ in LAYERS] + GROUPS:
        pk, fr = [], []
        for d in D:
            f, p = d["rate_spec"][name]
            pk.append(peak_freq(f, p, 5.0, 30.0))
            fr.append(band_fractions(f, p)["sigma 9-16"])
        lines.append(f"{LABEL.get(name, name):>10} | peak {_ms(pk, 1):>11} Hz | sigma fraction {_ms(fr, 3)}")
    return lines


def figures(D, out_prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1: spectra of the cumulative and layer traces, and of the groups
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for name, col in [("cortex", "k"), ("thalamus", "tab:red"), ("network", "tab:gray")] + \
            list(zip([n for n, _ in LAYERS], ["#1b4965", "#2471a3", "#5fa8d3", "#7f8c8d"])):
        f = D[0]["spec"][name][0]
        p = np.mean([d["spec"][name][1] for d in D], axis=0)
        ax[0].loglog(f[1:], p[1:], color=col, lw=2.2 if name in ("cortex", "thalamus", "network") else 1.2, label=name)
    ax[0].set_xlim(0.1, 40)
    ax[0].set_title("LFP proxy spectra: cumulative and layers (mean of runs)")
    ax[0].set_xlabel("frequency (Hz)")
    ax[0].set_ylabel("power (z²/Hz)")
    ax[0].legend(fontsize=8)
    cm = plt.get_cmap("tab10")
    for i, k in enumerate(GROUPS):
        f = D[0]["spec"][k][0]
        p = np.mean([d["spec"][k][1] for d in D], axis=0)
        ax[1].loglog(f[1:], p[1:], color=cm(i), lw=1.2, label=LABEL[k])
    ax[1].set_xlim(0.1, 40)
    ax[1].set_title("LFP proxy spectra by neuron group")
    ax[1].set_xlabel("frequency (Hz)")
    ax[1].legend(fontsize=7, ncol=2)
    for a in ax:
        a.axvspan(9, 16, color="gold", alpha=0.15)
        a.axvspan(0.3, 1.5, color="tab:blue", alpha=0.08)
    fig.tight_layout()
    fig.savefig(out_prefix + "_spectra.png", dpi=130)
    plt.close(fig)

    # 2: rates and UP-onset profiles
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    x = np.arange(len(GROUPS))
    for j, (key, col, lab) in enumerate([("mean", "tab:gray", "mean"), ("down", "tab:blue", "DOWN"), ("up", "tab:red", "UP")]):
        v = np.array([[d["groups"][k][key] for d in D] for k in GROUPS])
        ax[0].bar(x + (j - 1) * 0.27, v.mean(1), 0.27, yerr=v.std(1), color=col, label=lab, capsize=2)
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([LABEL[k] for k in GROUPS], rotation=45)
    ax[0].set_ylabel("rate (Hz/cell)")
    ax[0].set_yscale("log")
    ax[0].set_title("Firing rate by neuron group (mean ± SD over runs)")
    ax[0].legend()
    pt = D[0]["prof_t"]
    for i, k in enumerate(GROUPS):
        pr = np.mean([d["groups"][k]["profile"] for d in D], axis=0)
        pr = (pr - pr[(pt >= -400) & (pt < -100)].mean())
        pr = pr / (np.abs(pr).max() + 1e-12)
        ax[1].plot(pt, pr, color=cm(i), lw=1.4, label=LABEL[k])
    ax[1].axvline(0, color="k", lw=0.6)
    ax[1].set_xlim(-300, 500)
    ax[1].set_xlabel("time from UP onset (ms)")
    ax[1].set_ylabel("rate above DOWN level (normalised)")
    ax[1].set_title("UP-onset-locked rate by group")
    ax[1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_prefix + "_rates.png", dpi=130)
    plt.close(fig)

    # 3: example traces
    d = D[0]["trace"]
    names = ["network", "cortex", "thalamus"] + [n for n, _ in LAYERS] + GROUPS
    fig, ax = plt.subplots(figsize=(13, 9))
    tt = np.arange(len(d["cortex"])) / 1000.0 + 20.0
    for i, name in enumerate(names):
        y = d[name]
        ax.plot(tt, y / 4.0 - i, color="k" if name in ("network", "cortex", "thalamus") else "tab:blue", lw=0.7)
    ax.set_yticks([-i for i in range(len(names))])
    ax.set_yticklabels([LABEL.get(n, n) for n in names])
    ax.set_xlabel("time (s), seed 0")
    ax.set_title("LFP proxies (z-scored; cumulative traces in black)")
    fig.tight_layout()
    fig.savefig(out_prefix + "_traces.png", dpi=130)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True)
    ap.add_argument("--out-prefix", default="")
    ap.add_argument("--out-txt", default="")
    a = ap.parse_args(argv)
    D = collect(a.run)
    txt = "\n".join(tables(D))
    print(txt)
    if a.out_txt:
        open(a.out_txt, "w").write(txt + "\n")
    if a.out_prefix:
        figures(D, a.out_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
