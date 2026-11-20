"""
Membrane potential of the model cells: do they show UP and DOWN states?

Needs a run made with --vm-record N (ctx_thalamus_mpi.py): the npz then holds, for every
population, the mean soma potential over ALL its cells (vm_mean_<pop>, 1 kHz) and the
soma traces of the first N cells (vm_cells_<pop>). UP states are detected from the spikes
as elsewhere (stage2_report.up_onsets, L5E + L6E rate peaks).

Outputs (prefix given with --out-prefix):
  _traces.png  soma Vm of sampled L2/3 pyramidal cells, mean Vm of the pyramidal
               populations, L5E+L6E rate with the detected UP onsets
  _hist.png    Vm histograms (sampled L2/3 E cells; pooled), UP vs DOWN windows
  _stack.png   surface EEG proxy / deep EEG proxy / pyramidal-cell Vm, 7 s (classic layout)
  _onset.png   UP-onset-locked mean Vm of the populations (+-1.5 s) and of the sampled cells
plus a text table (--out-txt).

    python3 neuron/vm_report.py --run vm_s0.npz --out-prefix figures/vm --out-txt res/vm_table.txt
"""

import argparse

import numpy as np

import stage2_report as S

SKIP = 2000.0
UP_WIN, DOWN_WIN = (0.0, 300.0), (-400.0, -100.0)
POPS_E = ["l23e", "l4e", "l5e", "l6e"]


def load(path):
    z = np.load(path, allow_pickle=True)
    T = float(z["tstop"])
    ups = S.up_onsets(z["times"], z["gids"], z["ranges"].item(), T, SKIP)
    return z, T, np.asarray(ups, float)


def window_mean(v, ups, win, dt=1.0):
    out = []
    for u in ups:
        a, b = int((u + win[0]) / dt), int((u + win[1]) / dt)
        if a >= SKIP / dt and b < len(v):
            out.append(v[a:b].mean())
    return np.array(out)


def table(z, T, ups):
    dt = float(z["vm_t"][1] - z["vm_t"][0])
    out = ["Membrane potential, run of %.0f s, %d UP onsets (%.2f per s)"
           % (T / 1000, len(ups), len(ups) / ((T - SKIP) / 1000.0)), ""]
    out.append("%-6s %9s %9s %9s %12s" % ("pop", "mean Vm", "UP win", "DOWN win", "UP - DOWN (mV)"))
    for k in ["l23e", "l23i", "l4e", "l4i", "l5e", "l5i", "l6e", "l6i", "tc", "re"]:
        v = z["vm_mean_" + k].astype(float)
        up, dn = window_mean(v, ups, UP_WIN, dt), window_mean(v, ups, DOWN_WIN, dt)
        out.append("%-6s %9.1f %9.1f %9.1f %12.1f" % (k, v[int(SKIP / dt):].mean(), up.mean(),
                                                      dn.mean(), up.mean() - dn.mean()))
    out.append("")
    out.append("Sampled cells (soma Vm), UP window minus DOWN window (mV), mean over onsets:")
    for k in POPS_E:
        c = z["vm_cells_" + k].astype(float)
        d = [window_mean(c[i], ups, UP_WIN, dt).mean() - window_mean(c[i], ups, DOWN_WIN, dt).mean()
             for i in range(len(c))]
        lo = [np.percentile(c[i][int(SKIP / dt):], 5) for i in range(len(c))]
        hi = [np.percentile(c[i][int(SKIP / dt):], 95) for i in range(len(c))]
        out.append("  %-5s %s   5th-95th percentile of Vm: %.0f to %.0f mV (mean over cells)"
                   % (k, " ".join("%5.1f" % x for x in d), np.mean(lo), np.mean(hi)))
    return out


def figures(z, T, ups, pre):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = z["vm_t"]
    dt = float(t[1] - t[0])
    t0 = int(min(20, max(2, T / 1000.0 - 14))) * 1000
    w = (t >= t0) & (t < t0 + 12000)
    tt = t[w] / 1000.0

    # rate of L5E + L6E (UP detector input) in 5 ms bins
    R = z["ranges"].item()
    m = np.zeros(len(z["times"]), bool)
    n = 0
    for k in ("l5e", "l6e"):
        lo, hi = R[k]
        m |= (z["gids"] >= lo) & (z["gids"] < hi)
        n += hi - lo
    edges = np.arange(0, T + 5, 5.0)
    rate = np.histogram(z["times"][m], edges)[0] / n / 0.005
    rt = edges[:-1] / 1000.0

    c = z["vm_cells_l23e"].astype(float)
    nshow = min(5, len(c))
    fig, ax = plt.subplots(nshow + 2, 1, figsize=(14, 2 * (nshow + 2)), sharex=True,
                           gridspec_kw={"height_ratios": [1] * nshow + [1.3, 1]})
    for i in range(nshow):
        ax[i].plot(tt, c[i][w], "k", lw=0.6)
        ax[i].set_ylabel("L2/3 E\ncell %d\nVm (mV)" % i, fontsize=8)
    for k, col in (("l23e", "C0"), ("l4e", "C1"), ("l5e", "C3"), ("l6e", "C2")):
        ax[nshow].plot(tt, z["vm_mean_" + k][w], col, lw=0.9, label=k)
    ax[nshow].legend(ncol=4, frameon=False, fontsize=8, loc="upper right")
    ax[nshow].set_ylabel("mean Vm over\nall cells (mV)", fontsize=8)
    s = (rt >= t0 / 1000) & (rt < t0 / 1000 + 12)
    ax[nshow + 1].plot(rt[s], rate[s], "k", lw=0.6)
    ax[nshow + 1].set_ylabel("L5E+L6E\nrate (Hz)", fontsize=8)
    ax[nshow + 1].set_xlabel("time (s)")
    for a in ax:
        for u in ups:
            if t0 <= u < t0 + 12000:
                a.axvline(u / 1000.0, color="C3", lw=0.6, alpha=0.6)
    ax[0].set_title("Soma membrane potential of L2/3 pyramidal cells; red lines = detected UP onsets")
    fig.tight_layout()
    fig.savefig(pre + "_traces.png", dpi=120)
    plt.close(fig)

    # histograms
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    sk = int(SKIP / dt)
    for i in range(min(5, len(c))):
        ax[0].hist(c[i][sk:], bins=np.arange(-100, 50, 1), histtype="step", density=True,
                   label="cell %d" % i)
    ax[0].set_xlabel("soma Vm (mV)")
    ax[0].set_ylabel("density")
    ax[0].set_title("L2/3 E cells, whole run")
    ax[0].legend(frameon=False, fontsize=7)
    v = z["vm_mean_l23e"].astype(float)
    ax[1].hist(v[sk:], bins=60, color="C0")
    ax[1].set_xlabel("mean Vm of all L2/3 E cells (mV)")
    ax[1].set_title("population mean, whole run")
    up = np.concatenate([c[i][int(u / dt) + int(UP_WIN[0] / dt):int(u / dt) + int(UP_WIN[1] / dt)]
                         for i in range(len(c)) for u in ups
                         if SKIP < u < T - 400])
    dn = np.concatenate([c[i][int(u / dt) + int(DOWN_WIN[0] / dt):int(u / dt) + int(DOWN_WIN[1] / dt)]
                         for i in range(len(c)) for u in ups if SKIP < u < T - 400])
    ax[2].hist(dn, bins=np.arange(-100, 50, 1), alpha=0.6, density=True, label="DOWN window")
    ax[2].hist(up, bins=np.arange(-100, 50, 1), alpha=0.6, density=True, label="UP window")
    ax[2].set_xlabel("soma Vm (mV)")
    ax[2].legend(frameon=False)
    ax[2].set_title("L2/3 E cells: UP (0-300 ms) vs DOWN (-400...-100 ms) windows")
    fig.tight_layout()
    fig.savefig(pre + "_hist.png", dpi=120)
    plt.close(fig)

    # onset-locked
    h = int(1500 / dt)
    tt2 = np.arange(-h, h + 1) * dt / 1000.0
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for k, col in (("l23e", "C0"), ("l4e", "C1"), ("l5e", "C3"), ("l6e", "C2"), ("tc", "C4"),
                   ("re", "C5")):
        v = z["vm_mean_" + k].astype(float)
        A = np.array([v[int(u / dt) - h:int(u / dt) + h + 1] for u in ups
                      if u > SKIP + 1500 and u < T - 1600])
        ax[0].plot(tt2, A.mean(0) - A.mean(0)[: h // 3].mean(), col, label=k)
    ax[0].axvline(0, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel("time from UP onset (s)")
    ax[0].set_ylabel("mean Vm minus its value at -1.5...-1 s (mV)")
    ax[0].legend(frameon=False, ncol=2, fontsize=8)
    ax[0].set_title("population mean Vm locked to UP onsets")
    for i in range(min(5, len(c))):
        A = np.array([c[i][int(u / dt) - h:int(u / dt) + h + 1] for u in ups
                      if u > SKIP + 1500 and u < T - 1600])
        ax[1].plot(tt2, A.mean(0), lw=1, label="cell %d" % i)
    ax[1].axvline(0, color="k", ls=":", lw=0.8)
    ax[1].set_xlabel("time from UP onset (s)")
    ax[1].set_ylabel("soma Vm (mV)")
    ax[1].set_title("L2/3 E cells, onset-locked mean of single cells")
    ax[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(pre + "_onset.png", dpi=120)
    plt.close(fig)


def stack_figure(z, T, ups, pre, t0=None, span=7000.0, cell=0):
    """Layout of the classic three-trace figure (surface EEG, deep EEG, pyramidal cell).
    Surface EEG proxy: mean Vm of the superficial pyramidal cells (L2/3 E, L4 E), a
    depolarisation = surface-positive. Deep EEG proxy: mean Vm of the deep pyramidal
    cells (L5 E, L6 E) with the sign reversed, a depolarisation = depth-negative
    (Chauvette et al. 2010: the dipole polarity reverses between superficial and deep
    layers). Units: mV of mean membrane potential, NOT microvolts."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = z["vm_t"]
    if t0 is None:                       # window with the most UP onsets
        best = max(range(int(SKIP), int(T - span), 500),
                   key=lambda a: np.sum((ups >= a) & (ups < a + span)))
        t0 = best
    w = (t >= t0) & (t < t0 + span)
    tt = (t[w] - t0) / 1000.0

    def mean_of(keys):
        return np.mean([z["vm_mean_" + k].astype(float) for k in keys], axis=0)
    surf = mean_of(("l23e", "l4e"))
    deep = -mean_of(("l5e", "l6e"))
    fig, ax = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
    for a, y, lab in ((ax[0], surf, "surface EEG proxy\n(mean Vm L2/3+L4 E, mV)"),
                      (ax[1], deep, "deep EEG proxy\n(-mean Vm L5+L6 E, mV)"),
                      (ax[2], z["vm_cells_l23e"][cell].astype(float), "L2/3 pyramidal\ncell soma Vm (mV)"),
                      (ax[3], z["vm_cells_l5e"][cell].astype(float), "L5 pyramidal\ncell soma Vm (mV)")):
        a.plot(tt, y[w], "k", lw=0.8)
        a.set_ylabel(lab, fontsize=9)
        for u in ups:
            if t0 <= u < t0 + span:
                a.axvline((u - t0) / 1000.0, color="C3", lw=0.6, alpha=0.5)
    ax[3].set_xlabel("time (s); red = detected UP onsets")
    ax[0].set_title("Model: surface EEG proxy, deep EEG proxy and pyramidal-cell Vm (cf. classic "
                    "three-trace figure)")
    fig.tight_layout()
    fig.savefig(pre + "_stack.png", dpi=130)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--out-txt")
    a = ap.parse_args(argv)
    z, T, ups = load(a.run)
    txt = table(z, T, ups)
    print("\n".join(txt))
    if a.out_txt:
        open(a.out_txt, "w").write("\n".join(txt) + "\n")
    figures(z, T, ups, a.out_prefix)
    stack_figure(z, T, ups, a.out_prefix)


if __name__ == "__main__":
    main()
