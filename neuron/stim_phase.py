"""
Offline SO-phase and spindle analysis for the stimulation experiments
(PLAN-stimulation-experiments.md, WP2).

Per tone (stim arm: real tone; sham arm: fake tone) and per arm:
  * SO phase at the tone, 4 bins: rising / peak / falling / trough
  * spindle present at the tone, and its phase relative to the SO
  * responses in the 1 s after the tone: SO amplitude, P(spindle), spindle
    power, SO phase of the spindle peak

Conventions
  SO signal    cortical E rate (L2/3, L4, L5, L6), 10 ms bins, band 0.3-1.5 Hz.
  SO phase     Hilbert phase of that band: 0 = maximum of the UP-state
               activity ("peak"), -90 deg = rising (DOWN -> UP), +90 deg =
               falling, 180 deg = trough (DOWN).
  spindle      RE volley train (volley_stats), >= 6 volleys; present at t if the
               last volley <= t is within TRAIN_MAX_IVI and the train so far has
               >= 3 volleys. Spindle phase = SO phase at the train's RE
               participation maximum *so far* (the final maximum may lie after
               the tone, and in the stim arm that is the tone's own effect).

Why a causal phase estimate: the arms diverge after the first tone, so a tone
in the stim arm cannot be classified from the sham run, and an acausal
(filtfilt/Hilbert) phase at the tone leaks the tone's own response. Both arms
are therefore classified with `phase_causal`, which uses data before the tone
only (symmetric FIR + AR forecast). `--validate` measures its error against
the acausal phase on unperturbed (sham) runs. Tuned on 3 sham runs: longer
FIRs (1.5 s half-length) were worse (median error 56 deg, 41% same bin) than
this one (25 deg, ~60% same bin); it is a 0.3-1.5 Hz low-pass in effect.

    python3 neuron/stim_phase.py --validate res/.../cl_sham_s0.npz ...
    python3 neuron/stim_phase.py --pair stim.npz=sham.npz [--pair ...] \\
        --out-png figures/stim_graph_a.png --rows rows.pkl
"""

import argparse
import pickle
import sys

import numpy as np
from scipy.signal import butter, firwin, hilbert, sosfiltfilt

import cl_analyze as CLA
import volley_stats as V

BIN_MS = 10.0
FS = 1000.0 / BIN_MS
SO_BAND = (0.3, 1.5)
SP_BAND = (10.0, 15.0)
CX = ("l23e", "l4e", "l5e", "l6e")
SO_NAMES = ("rising", "peak", "falling", "trough")
SP_NAMES = ("none",) + SO_NAMES
WIN_MS = 1000.0


# ---------------------------------------------------------------- SO signal
def cortex_rate(npz):
    """Cortical E population rate, Hz/cell, 10 ms bins (bin k = [k, k+1)*BIN_MS)."""
    t, g, R, T = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    m = np.zeros_like(t, bool)
    n = 0
    for k in CX:
        lo, hi = R[k]
        m |= (g >= lo) & (g < hi)
        n += hi - lo
    bins = np.arange(0, T + BIN_MS, BIN_MS)
    return np.histogram(t[m], bins)[0] / n / (BIN_MS / 1000.0)


def so_filter(x):
    return sosfiltfilt(butter(2, SO_BAND, "bp", fs=FS, output="sos"), x - x.mean())


def phase_acausal(x):
    """(phase, amplitude) of the SO band at every bin; zero-phase, uses the
    whole run, so valid only where the run is unperturbed."""
    a = hilbert(so_filter(x))
    return np.angle(a), np.abs(a)


_TAPS = firwin(2 * 25 + 1, SO_BAND, pass_zero=False, fs=FS)


def _ar_forecast(y, n, p=10):
    X = np.column_stack([y[i:len(y) - p + i] for i in range(p)])
    c = np.linalg.lstsq(X, y[p:], rcond=None)[0]
    buf = list(y[-p:])
    for _ in range(n):
        buf.append(float(np.dot(c, buf[-p:])))
    return np.array(buf[p:])


def phase_causal(x, t_ms, win_s=30.0, ds=5, p=10):
    """SO phase at time t_ms from bins that end at or before t_ms only.
    Symmetric FIR band-pass (0.5 s half-length) on the last win_s seconds,
    valid up to t - 0.25 s; the last 0.25 s are forecast with an AR(p) model fitted on the filtered
    (decimated) series; phase from the Hilbert transform of the extended series.
    Returns nan if there is less than ~8 s of history."""
    i = int(t_ms // BIN_MS)
    h = (len(_TAPS) - 1) // 2
    seg = x[max(0, i - int(win_s * FS)):i]
    if len(seg) < int(8 * FS):
        return np.nan
    y = np.convolve(seg - seg.mean(), _TAPS, "valid")
    yd = y[::-1][::ds][::-1]
    lag = int(round((h + 1) / ds))
    fc = _ar_forecast(yd, lag + 15, p)
    a = hilbert(np.concatenate([yd, fc]))
    return float(np.angle(a[len(yd) - 1 + lag]))


def so_bin(ph):
    """0 rising [-135,-45), 1 peak [-45,45), 2 falling [45,135), 3 trough."""
    if not np.isfinite(ph):
        return -1
    return int(((ph + 0.75 * np.pi) % (2 * np.pi)) // (0.5 * np.pi))


# ------------------------------------------------------------------ spindles
def volleys_of(npz, skip=0.0):
    t, g, R, T = npz["times"], npz["gids"], npz["ranges"].item(), float(npz["tstop"])
    return V.re_volleys(t, g, R, T, skip_ms=skip)


def trains_of(volleys, min_cycles=6):
    """Spindle trains with their participation maximum (t_peak)."""
    out, cur = [], []

    def close():
        if len(cur) >= min_cycles:
            k = int(np.argmax([v["re_part"] for v in cur]))
            out.append({"t0": cur[0]["t"], "t1": cur[-1]["t"], "cycles": len(cur),
                        "t_peak": cur[k]["t"]})
    for v in volleys:
        if cur and v["t"] - cur[-1]["t"] > V.TRAIN_MAX_IVI:
            close()
            cur = []
        cur.append(v)
    if cur:
        close()
    return out


def spindle_state(volleys, t_ms, min_so_far=3):
    """Causal: (present, t_peak_so_far) from the volleys at or before t_ms."""
    vt = np.array([v["t"] for v in volleys])
    k = int(np.searchsorted(vt, t_ms, "right")) - 1
    if k < 0 or t_ms - vt[k] > V.TRAIN_MAX_IVI:
        return False, np.nan
    j = k
    while j > 0 and vt[j] - vt[j - 1] <= V.TRAIN_MAX_IVI:
        j -= 1
    if k - j + 1 < min_so_far:
        return False, np.nan
    part = [volleys[q]["re_part"] for q in range(j, k + 1)]
    return True, float(vt[j + int(np.argmax(part))])


# ------------------------------------------------------------------- a run
class Run:
    """Everything derived once per npz."""

    def __init__(self, path):
        self.z = np.load(path, allow_pickle=True)
        self.x = cortex_rate(self.z)
        self.phase_a, self.amp_a = phase_acausal(self.x)
        self.xf = so_filter(self.x)
        self.volleys = volleys_of(self.z)
        self.trains = trains_of(self.volleys)
        R, T = self.z["ranges"].item(), float(self.z["tstop"])
        lo, hi = R["re"]
        tt, gg = self.z["times"], self.z["gids"]
        m = (gg >= lo) & (gg < hi)
        n1 = int(T) + 1
        c = np.histogram(tt[m], np.arange(0, n1 + 1.0, 1.0))[0] / (hi - lo) * 1000.0
        a = hilbert(sosfiltfilt(butter(4, SP_BAND, "bp", fs=1000.0, output="sos"), c))
        self.sp_env2 = np.abs(a) ** 2
        self.T = T

    def at(self, t_ms):
        return min(int(t_ms // BIN_MS), len(self.phase_a) - 1)

    def classify(self, t_ms):
        ph = phase_causal(self.x, t_ms)
        pres, tp = spindle_state(self.volleys, t_ms)
        sp_ph = float(self.phase_a[self.at(tp)]) if pres else np.nan  # SO phase at the peak
        return {"t": float(t_ms), "so_phase": ph, "so_bin": so_bin(ph),
                "so_phase_acausal": float(self.phase_a[self.at(t_ms)]),
                "sp_present": bool(pres), "sp_phase": sp_ph,
                "sp_bin": (so_bin(sp_ph) + 1) if pres else 0}

    def respond(self, t_ms, w=WIN_MS):
        a, b = self.at(t_ms), self.at(t_ms + w)
        peaks = [tr for tr in self.trains if t_ms < tr["t_peak"] <= t_ms + w]
        return {"so_amp": float(np.sqrt(np.mean(self.xf[a:b] ** 2))),
                "sp_any": bool(peaks),
                "sp_pow": float(np.mean(self.sp_env2[int(t_ms):int(t_ms + w)])),
                "sp_phase_after": float(self.phase_a[self.at(peaks[0]["t_peak"])]) if peaks else np.nan}


def tone_rows(run, tones, skip=2000.0, w=WIN_MS):
    rows = []
    for t in tones:
        if t < skip or t > run.T - w - 100.0:
            continue
        r = run.classify(t)
        r.update(run.respond(t, w))
        rows.append(r)
    return rows


# -------------------------------------------------------------- validation
def validate(paths, n=400, seed=0):
    rng = np.random.default_rng(seed)
    conf = np.zeros((4, 4), int)
    print(f"{'run':<22}{'UP onsets: phase mean / R':>28}{'causal err (deg): med / p90':>30}"
          f"{'same bin':>10}{'spindle time':>14}{'sp-peak phase mean / R':>26}")
    for p in paths:
        r = Run(p)
        import stage2_report as S
        z = r.z
        ups = S.up_onsets(z["times"], z["gids"], z["ranges"].item(), r.T, 2000.0)
        pu = r.phase_a[[r.at(u) for u in ups]]
        mu = np.angle(np.mean(np.exp(1j * pu)))
        Ru = np.abs(np.mean(np.exp(1j * pu)))
        ts = rng.uniform(10000.0, r.T - 2000.0, n)
        err, same = [], []
        for t in ts:
            pc, pa = phase_causal(r.x, t), r.phase_a[r.at(t)]
            err.append(abs(np.angle(np.exp(1j * (pc - pa)))))
            same.append(so_bin(pc) == so_bin(pa))
            conf[so_bin(pa), so_bin(pc)] += 1
        err = np.degrees(err)
        occ = np.mean([spindle_state(r.volleys, t)[0] for t in ts])
        sp = np.array([r.phase_a[r.at(tr["t_peak"])] for tr in r.trains])
        ms, Rs = (np.degrees(np.angle(np.mean(np.exp(1j * sp)))), np.abs(np.mean(np.exp(1j * sp)))) \
            if len(sp) else (np.nan, np.nan)
        print(f"{p.split('/')[-1]:<22}{np.degrees(mu):>14.0f} / {Ru:4.2f}{np.median(err):>16.0f} / {np.percentile(err, 90):4.0f}"
              f"{np.mean(same):>10.0%}{occ:>14.0%}{ms:>18.0f} / {Rs:4.2f}")
    print("\nbin confusion, rows = true (acausal), cols = causal estimate; order "
          + ", ".join(SO_NAMES) + f"\n{conf}")


# ----------------------------------------------------------------- Graph A
def collect_pairs(pairs, skip=2000.0):
    stim, sham = [], []
    for sp, sh in pairs:
        for path, rows in ((sp, stim), (sh, sham)):
            run = Run(path)
            t, _ = CLA.tones_of(run.z)
            rows.extend(tone_rows(run, t, skip))
    return stim, sham


def grid(rows, key, fn=np.nanmean):
    """(4 SO bins) x (5 spindle bins) array of fn(rows[key]) and counts."""
    G = np.full((4, 5), np.nan)
    N = np.zeros((4, 5), int)
    for i in range(4):
        for j in range(5):
            v = [r[key] for r in rows if r["so_bin"] == i and r["sp_bin"] == j]
            N[i, j] = len(v)
            if len(v):
                G[i, j] = fn(np.array(v, float))
    return G, N


def graph_a(stim, sham, out_png, min_n=5):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    meas = [("so_amp", "SO amplitude 0-1 s", np.nanmean),
            ("sp_any", "P(spindle peak in 0-1 s)", np.nanmean),
            ("sp_pow", "spindle power 0-1 s (RE 10-15 Hz)", np.nanmean)]
    fig, axs = plt.subplots(len(meas), 3, figsize=(12, 3.0 * len(meas)))
    for row, (key, title, fn) in zip(axs, meas):
        A, Na = grid(stim, key, fn)
        B, Nb = grid(sham, key, fn)
        D = A - B
        bad = (Na < min_n) | (Nb < min_n)
        for ax, M, lab in ((row[0], np.where(Na >= min_n, A, np.nan), "stim"),
                           (row[1], np.where(Nb >= min_n, B, np.nan), "sham"),
                           (row[2], np.where(bad, np.nan, D), "stim - sham")):
            im = ax.imshow(M, cmap="viridis" if lab != "stim - sham" else "RdBu_r",
                           aspect="auto", origin="lower")
            ax.set_title(f"{title}: {lab}", fontsize=8)
            ax.set_xticks(range(5))
            ax.set_xticklabels(SP_NAMES, fontsize=7)
            ax.set_yticks(range(4))
            ax.set_yticklabels(SO_NAMES, fontsize=7)
            plt.colorbar(im, ax=ax, fraction=0.046)
            for i in range(4):
                for j in range(5):
                    n = (Na[i, j], Nb[i, j]) if lab == "stim - sham" else ((Na if lab == "stim" else Nb)[i, j],)
                    ax.text(j, i, "/".join(str(k) for k in n), ha="center", va="center",
                            fontsize=6, color="w" if lab != "stim - sham" else "k")
    axs[-1][0].set_xlabel("spindle phase relative to SO (none = no spindle)", fontsize=8)
    for row in axs:
        row[0].set_ylabel("SO phase at tone", fontsize=8)
    fig.suptitle(f"Graph A (grid view): cells with n < {min_n} in either arm are blank; "
                 "numbers = tones per cell", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validate", nargs="+", metavar="NPZ")
    ap.add_argument("--pair", action="append", help="stim.npz=sham.npz (closed-loop run pair)")
    ap.add_argument("--rows", default="", help="pickle for the per-tone rows")
    ap.add_argument("--out-png", default="")
    ap.add_argument("--min-n", type=int, default=5)
    a = ap.parse_args(argv)
    if a.validate:
        validate(a.validate)
    if a.pair:
        stim, sham = collect_pairs([p.split("=", 1) for p in a.pair])
        if a.rows:
            pickle.dump((stim, sham), open(a.rows, "wb"))
        print(f"{len(stim)} stim tones, {len(sham)} sham tones")
        for lab, rows in (("stim", stim), ("sham", sham)):
            _, N = grid(rows, "so_amp")
            print(f"{lab} tones per cell (rows {SO_NAMES}, cols {SP_NAMES}):\n{N}")
        if a.out_png:
            graph_a(stim, sham, a.out_png, a.min_n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
