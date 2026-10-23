"""
Frequency comparison of stimulus trains (PLAN-stimulation-experiments.md,
"Stimulation rates from the literature"): stim vs sham for trains of n stimuli
at rate 1000/isi Hz, from --open-loop runs with --ol-n and --ol-isi lists.

Per protocol (n, isi) and per pulse index k (response fading along the train):
  TC spikes/cell in 50 ms, P(UP onset within 300 ms), P(spindle within 300 ms).
Per train, over its span [t0, t_last + isi]:
  UP-onset rate, and entrainment: the vector strength R of the UP onsets' phase
  in the stimulation period (phase = 2 pi ((t_up - t0) mod isi) / isi), pooled
  over trains, Rayleigh p.  Sham fake pulses at the same times are the control.
Per train, after it ([t_last + 1 s, t_last + 6 s]): UP rate and spindle rate,
the longer-term effect that the literature reports as absent.

    python3 neuron/stim_train.py --pair stim.npz=sham.npz [...] --cache rows.pkl \\
        --out-prefix figures/f1 --out-txt table.txt
"""

import argparse
import os
import pickle

import numpy as np

import aud_phase as AP
import stage2_report as S
import stim_phase as SPH
import volley_stats as V


def pulses_of(z):
    """list of (t, protocol (n, isi), train index, pulse index k)."""
    out = []
    for j, (t0, n, isi) in enumerate(zip(z["ol_t"], z["ol_n"], z["ol_isi"])):
        for k in range(int(n)):
            out.append((float(t0 + k * isi), (int(n), float(isi)), j, k))
    return out


def run_rows(path, skip=2000.0, post=(1000.0, 6000.0)):
    """(pulse rows, train rows) of one run."""
    z = np.load(path, allow_pickle=True)
    T = float(z["tstop"])
    P = pulses_of(z)
    meas = {round(r["t0"], 3): r for r in AP.measure(z, [p[0] for p in P], skip, 2.0)}
    pulse_rows = []
    for t, proto, j, k in P:
        r = meas.get(round(t, 3))
        if r is None:
            continue
        pulse_rows.append({"t": t, "n": proto[0], "isi": proto[1], "train": j, "k": k,
                           "tc": r["tc"], "evoked_up": r["evoked_up"], "evoked_sp": r["evoked_sp"]})
    t_, g_, R = z["times"], z["gids"], z["ranges"].item()
    ups = S.up_onsets(t_, g_, R, T, skip)
    trains = V.volley_trains(V.re_volleys(t_, g_, R, T, skip_ms=skip))
    sp_t0 = np.array([tr["t0"] for tr in trains if tr["cycles"] >= 6])
    train_rows = []
    for j, (t0, n, isi) in enumerate(zip(z["ol_t"], z["ol_n"], z["ol_isi"])):
        end = t0 + (n - 1) * isi
        a, b = t0, end + isi
        if t0 < skip or b > T or end + post[1] > T:
            continue
        u = ups[(ups >= a) & (ups < b)]
        ph = 2 * np.pi * (((u - t0) % isi) / isi)
        pu = ups[(ups >= end + post[0]) & (ups < end + post[1])]
        ps = sp_t0[(sp_t0 >= end + post[0]) & (sp_t0 < end + post[1])]
        train_rows.append({"train": j, "n": int(n), "isi": float(isi), "span": float(b - a),
                           "n_up": len(u), "phases": ph,
                           "post_up_hz": len(pu) / ((post[1] - post[0]) / 1000.0),
                           "post_sp_hz": len(ps) / ((post[1] - post[0]) / 1000.0)})
    return pulse_rows, train_rows


def collect(pairs):
    out = {"stim": ([], []), "sham": ([], [])}
    for k, pair in enumerate(pairs):
        for arm, path in zip(("stim", "sham"), pair):
            pr, tr = run_rows(path)
            for r in pr + tr:
                r["seed"] = k
            out[arm][0].extend(pr)
            out[arm][1].extend(tr)
    return out


def rayleigh(phases):
    return SPH.circ_summary(phases)


def protocols(D):
    return sorted({(r["n"], r["isi"]) for arm in D.values() for r in arm[1]})


def summary(D, min_n=1):
    lines = [f"{'protocol':>14} | {'trains s/h':>10} | {'UP Hz in train':>16} | {'UP vector strength R':>26} | "
             f"{'UP Hz after':>14} | {'spindle Hz after':>17}"]
    for n, isi in protocols(D):
        if n < 2:
            continue
        a = [r for r in D["stim"][1] if (r["n"], r["isi"]) == (n, isi)]
        b = [r for r in D["sham"][1] if (r["n"], r["isi"]) == (n, isi)]
        if len(a) < min_n or len(b) < min_n:
            continue
        up = lambda rs: sum(r["n_up"] for r in rs) / (sum(r["span"] for r in rs) / 1000.0)
        pa = np.concatenate([r["phases"] for r in a]) if a else np.empty(0)
        pb = np.concatenate([r["phases"] for r in b]) if b else np.empty(0)
        _, Ra, pva, na = rayleigh(pa)
        _, Rb, pvb, nb = rayleigh(pb)
        m = lambda rs, k: np.mean([r[k] for r in rs])
        lines.append(f"{n:>3} @ {1000/isi:4.2f} Hz | {len(a):>4}/{len(b):<5} | {up(a):5.2f} vs {up(b):4.2f}   | "
                     f"R {Ra:4.2f} (n={na}, p={pva:.3f}) vs {Rb:4.2f} (p={pvb:.3f}) | "
                     f"{m(a,'post_up_hz'):5.2f} vs {m(b,'post_up_hz'):4.2f} | {m(a,'post_sp_hz'):5.2f} vs {m(b,'post_sp_hz'):4.2f}")
    return lines


def figures(D, out_prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    isis = sorted({r["isi"] for r in D["stim"][0] if r["n"] >= 2})
    cmap = plt.get_cmap("viridis")
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, (key, title) in zip(axs, (("tc", "TC spikes/cell, 50 ms after the pulse"),
                                      ("evoked_up", "P(UP onset within 300 ms)"),
                                      ("evoked_sp", "P(spindle within 300 ms)"))):
        for i, isi in enumerate(isis):
            for arm, ls in (("stim", "-"), ("sham", ":")):
                ks = sorted({r["k"] for r in D[arm][0] if r["isi"] == isi})
                m = [np.mean([r[key] for r in D[arm][0] if r["isi"] == isi and r["k"] == k]) for k in ks]
                ax.plot(np.array(ks) + 1, m, ls, color=cmap(i / max(len(isis) - 1, 1)), marker="o", ms=3,
                        label=f"{1000/isi:.2f} Hz {arm}" if arm == "stim" or i == 0 else None)
        ax.set_xlabel("pulse index in the train")
        ax.set_title(title, fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    axs[0].legend(frameon=False, fontsize=6, ncol=2)
    fig.suptitle("Response along the train, by stimulation rate (solid = stim, dotted = sham fake pulses)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_prefix + "_per_pulse.png", dpi=150)

    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    prot = [p for p in protocols(D) if p[0] >= 2]
    labels = [f"{n}@{1000/i:.2f}" for n, i in prot]
    x = np.arange(len(prot))
    for ax, key in zip(axs, ("up", "R", "post")):
        for arm, col, off in (("sham", "tab:gray", -0.1), ("stim", "tab:red", 0.1)):
            ys = []
            for n, isi in prot:
                rs = [r for r in D[arm][1] if (r["n"], r["isi"]) == (n, isi)]
                if not rs:
                    ys.append(np.nan)
                elif key == "up":
                    ys.append(sum(r["n_up"] for r in rs) / (sum(r["span"] for r in rs) / 1000.0))
                elif key == "R":
                    ys.append(rayleigh(np.concatenate([r["phases"] for r in rs]))[1])
                else:
                    ys.append(np.mean([r["post_up_hz"] for r in rs]))
            ax.plot(x + off, ys, "o-", color=col, label=arm)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, fontsize=7)
        ax.set_title({"up": "UP-onset rate during the train (Hz)",
                      "R": "UP vector strength in the stimulus period",
                      "post": "UP rate 1-6 s after the train (Hz)"}[key], fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    axs[0].legend(frameon=False, fontsize=8)
    axs[1].set_xlabel("protocol: n stimuli @ rate (Hz)", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_prefix + "_entrainment.png", dpi=150)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", action="append", required=True)
    ap.add_argument("--cache", default="")
    ap.add_argument("--out-prefix", default="")
    ap.add_argument("--out-txt", default="")
    a = ap.parse_args(argv)
    if a.cache and os.path.exists(a.cache):
        D = pickle.load(open(a.cache, "rb"))
    else:
        D = collect([p.split("=", 1) for p in a.pair])
        if a.cache:
            pickle.dump(D, open(a.cache, "wb"))
    txt = "\n".join(summary(D))
    print(txt)
    if a.out_txt:
        open(a.out_txt, "w").write(txt + "\n")
    if a.out_prefix:
        figures(D, a.out_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
