"""
MPI-parallel NEURON thalamic network for MareNostrum 5.

`tc_network_nrn.ThalamicNet` is serial (one process, all cells). This module
rebuilds the same biophysics on NEURON's ParallelContext so the population can
be distributed over MPI ranks -- the prerequisite for running the HH + SK2 +
Ca-dependent model at column scale instead of the 20-cell laptop network.

Cells are the SAME objects as the serial model (tc_neuron.TCCell / RECell), so
the mechanism set is identical: hh2 (Traub-Miles HH), itd (Destexhe T-current),
cad (submembrane Ca pool), sk2 (Ca-activated K), ihca (Ca-dependent I_h).

Run:
    # scaling benchmark -- the measurement that decides the architecture
    srun -n 100 .venv-neuron/bin/nrniv -python -mpi tc_mpi.py --bench

    # production
    srun -n 100 .venv-neuron/bin/nrniv -python -mpi tc_mpi.py \
         --n-tc 2500 --n-re 2500 --tstop 200000 --out out/mn5_nrn.h5

LIMITATION -- gap junctions are rank-local in this version. NEURON needs
pc.setup_transfer()/source_var()/target_var() for cross-rank electrical
coupling; here RE cells are block-assigned so each rank holds a contiguous
group and gap junctions are made only within a rank. Electrical coupling is
therefore clustered rather than global. Flagged rather than silently dropped:
see neuron/README.md.
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from neuron import h                                    # noqa: E402

# MUST come before ParallelContext when launched as `mpirun -n N python ...`.
# Plain python does not initialise NEURON's MPI, so every process would report
# nhost=1, believe it is rank 0, and run the WHOLE network independently -- N
# identical serial copies silently consuming N cpus. nrniv -python -mpi does
# this for us; calling it explicitly makes both launch paths correct.
h.nrnmpi_init()

import tc_neuron as T                                    # noqa: E402

h.load_file("stdrun.hoc")


class ParallelThalamicNet:
    """The TC<->RE loop distributed over MPI ranks.

    gids 0..n_tc-1        -> TC (relay) cells
    gids n_tc..n_tc+n_re-1 -> RE (reticular) cells

    RE cells are BLOCK-assigned (contiguous per rank) so intra-rank gap
    junctions form local clusters; TC cells are round-robin for load balance.
    """

    def __init__(self, n_tc=10, n_re=10, seed=1,
                 g_re_tc=0.015, g_tc_re=0.004, g_re_re=0.002, g_gap=0.03,
                 g_cort=0.01, g_cort_tc=0.035, so_freq=1.0,
                 tc_e_pas=-80.0, re_e_pas=-82.0,
                 gsk_re=5e-5, gh_tc=4e-4, conv=10):
        self.pc = h.ParallelContext()
        self.rank = int(self.pc.id())
        self.nhost = int(self.pc.nhost())
        self.n_tc, self.n_re = n_tc, n_re
        self.so_freq = so_freq
        self.rng = np.random.default_rng(seed + self.rank)
        self.cells, self.syns, self.ncs, self.stims, self.gaps = {}, {}, [], [], []

        # ---- distribute gids -------------------------------------------------
        self.tc_gids = list(range(self.rank, n_tc, self.nhost))        # round-robin
        per = int(np.ceil(n_re / self.nhost))                          # block
        lo = n_tc + self.rank * per
        self.re_gids = list(range(lo, min(lo + per, n_tc + n_re)))

        for gid in self.tc_gids:
            c = T.TCCell(gsk=0.0, gh=gh_tc)
            c.soma.e_pas = tc_e_pas
            self._register(gid, c)
        for gid in self.re_gids:
            c = T.RECell(gsk=gsk_re)
            c.soma.e_pas = re_e_pas
            self._register(gid, c)

        # ---- synapses (created on the TARGET's rank) -------------------------
        # Synaptic kinetics COPIED from the serial model's _syn_connect calls
        # (tc_network_nrn.py) -- these must match or the two networks are not
        # the same model.
        for gid in self.tc_gids:                    # RE -> TC : GABA_A
            s = h.Exp2Syn(self.cells[gid].soma(0.5))
            s.tau1, s.tau2, s.e = 1.0, 8.0, -85.0
            self.syns[("gabaa", gid)] = s
            s2 = h.Exp2Syn(self.cells[gid].soma(0.5))   # cortex -> TC : AMPA
            s2.tau1, s2.tau2, s2.e = 0.5, 2.0, 0.0
            self.syns[("ampa_c", gid)] = s2
        for gid in self.re_gids:                    # TC -> RE : AMPA
            s = h.Exp2Syn(self.cells[gid].soma(0.5))
            s.tau1, s.tau2, s.e = 0.5, 2.0, 0.0
            self.syns[("ampa", gid)] = s
            s2 = h.Exp2Syn(self.cells[gid].soma(0.5))   # RE -> RE : GABA_A
            s2.tau1, s2.tau2, s2.e = 1.0, 6.0, -75.0
            self.syns[("gabaa_re", gid)] = s2
            s3 = h.Exp2Syn(self.cells[gid].soma(0.5))   # cortex -> RE : AMPA
            s3.tau1, s3.tau2, s3.e = 0.5, 2.0, 0.0
            self.syns[("ampa_c", gid)] = s3

        # ---- wire across ranks via gid_connect ------------------------------
        # Rule RECONCILED with the serial model (tc_network_nrn.py):
        #   RE -> TC : each TC draws n_re//2 reticular sources, weight g/k
        #   TC -> RE : each RE draws n_tc//2 relay sources,     weight g/k
        #   RE -> RE : nearest neighbours |i-j| <= 2,           weight g (NOT
        #              normalised, matching _wire_re_re)
        # The g/k FAN-IN NORMALISATION is the critical part: an earlier version
        # used a fixed conv=10 at full weight, giving each target ~10x the
        # intended total conductance and running the loop at 20 Hz instead of
        # 13 Hz.
        #
        # The serial model draws from one shared RNG stream, so its exact
        # realisation depends on call order and cannot be reproduced
        # rank-independently. Seeding per TARGET gid instead gives the same
        # convergence and weight distribution, reproducible and independent of
        # rank count -- statistically equivalent, not bit-identical.
        k_re = max(1, n_re // 2)
        k_tc = max(1, n_tc // 2)
        for gid in self.tc_gids:                                   # RE -> TC
            for src in self._draw(gid, n_tc, n_tc + n_re, k_re):
                self._connect(src, self.syns[("gabaa", gid)], g_re_tc / k_re, 1.0)
        for gid in self.re_gids:                                   # TC -> RE
            for src in self._draw(gid, 0, n_tc, k_tc):
                self._connect(src, self.syns[("ampa", gid)], g_tc_re / k_tc, 1.0)
            i = gid - n_tc                                         # RE -> RE
            for j in range(max(0, i - 2), min(n_re, i + 3)):
                if j != i:
                    self._connect(n_tc + j, self.syns[("gabaa_re", gid)],
                                  g_re_re, 1.0)

        # ---- rank-local gap junctions between RE cells -----------------------
        if g_gap > 0:
            loc = [self.cells[g] for g in self.re_gids]
            for a, b in zip(loc[:-1], loc[1:]):
                self.gaps.append(T.gap_junction(a, b, g=g_gap))

        # ---- corticothalamic drive (UP-state gated burst volleys) ------------
        self._drive(g_cort, g_cort_tc)

        # ---- spike recording -------------------------------------------------
        self.tspk, self.gspk = h.Vector(), h.Vector()
        self.pc.spike_record(-1, self.tspk, self.gspk)

    # ------------------------------------------------------------------ utils
    def _register(self, gid, cell):
        self.cells[gid] = cell
        self.pc.set_gid2node(gid, self.rank)
        nc = h.NetCon(cell.soma(0.5)._ref_v, None, sec=cell.soma)
        nc.threshold = -10
        self.pc.cell(gid, nc)
        self.ncs.append(nc)

    def _draw(self, seed_gid, lo, hi, k):
        r = np.random.default_rng(seed_gid)
        n = hi - lo
        return (lo + r.choice(n, size=min(k, n), replace=False)) if n else []

    def _connect(self, src_gid, syn, weight, delay):
        nc = self.pc.gid_connect(int(src_gid), syn)
        nc.weight[0], nc.delay = weight, delay
        self.ncs.append(nc)

    def _drive(self, g_cort, g_cort_tc, tstop_hint=200000.0):
        """One burst volley per slow-oscillation cycle, onto RE and TC."""
        period = 1000.0 / self.so_freq
        n_cyc = int(tstop_hint / period) + 1
        for k in range(n_cyc):
            t0 = 300.0 + k * period
            ns = h.NetStim()
            ns.interval, ns.number, ns.start, ns.noise = 4.0, 5, t0, 0.0
            self.stims.append(ns)
            for gid in self.re_gids:
                nc = h.NetCon(ns, self.syns[("ampa_c", gid)])
                nc.weight[0], nc.delay = g_cort, 1.0
                self.ncs.append(nc)
            for gid in self.tc_gids:
                nc = h.NetCon(ns, self.syns[("ampa_c", gid)])
                nc.weight[0], nc.delay = g_cort_tc, 1.0
                self.ncs.append(nc)

    # -------------------------------------------------------------------- run
    def run(self, tstop=12000.0, celsius=36.0, dt=0.025):
        self.pc.set_maxstep(10)
        h.celsius, h.dt = celsius, dt
        h.finitialize(-75)
        t0 = time.time()
        self.pc.psolve(tstop)
        wall = time.time() - t0
        wall = self.pc.allreduce(wall, 2)          # max over ranks
        return wall

    def gather(self):
        """Collect all spikes onto rank 0. Returns (times, gids) or (None, None)."""
        t_all = self.pc.py_gather(np.asarray(self.tspk), 0)
        g_all = self.pc.py_gather(np.asarray(self.gspk), 0)
        if self.rank != 0:
            return None, None
        t = np.concatenate(t_all) if t_all else np.array([])
        g = np.concatenate(g_all) if g_all else np.array([])
        o = np.argsort(t)
        return t[o], g[o]


def _report(net, t, g, wall, tstop):
    """Rank-0 summary: frequency measured with an UNCONSTRAINED spectral search."""
    n_tc = net.n_tc
    print("\n=== NEURON MPI thalamic network ===")
    print("ranks %d | TC %d | RE %d | tstop %.0f ms | wall %.1f s (%.2f x realtime)"
          % (net.nhost, n_tc, net.n_re, tstop, wall, wall / (tstop / 1000.0)))
    if t is None or len(t) < 20:
        print("too few spikes to analyse (%d)" % (0 if t is None else len(t)))
        return
    m = t > 2000.0                                   # discard startup transient
    t, g = t[m], g[m]
    print("spikes %d  (TC %d, RE %d)" % (len(t), int((g < n_tc).sum()),
                                         int((g >= n_tc).sum())))
    brk = np.where(np.diff(t) > 40.0)[0]
    st = np.concatenate(([t[0]], t[brk + 1]))
    if len(st) > 2:
        ibi = np.diff(st)
        med = float(np.median(ibi))
        iqr = float(np.percentile(ibi, 75) - np.percentile(ibi, 25))
        print("population burst rate %.2f Hz (IQR %.0f ms, n=%d)"
              % (1000.0 / med, iqr, len(st)))
        if abs(med - 1000.0 / net.so_freq) < 30.0 and iqr < 25.0:
            print("  WARNING: drive-locked -- this is the %.1f Hz drive, not an "
                  "intrinsic rhythm" % net.so_freq)
    # unconstrained spectrum (the discipline that caught the 12 Hz error)
    bins = np.arange(2000.0, t.max() + 1.0, 1.0)
    rate = np.histogram(t, bins=bins)[0].astype(float)
    rate -= rate.mean()
    f = np.fft.rfftfreq(len(rate), 1e-3)
    p = np.abs(np.fft.rfft(rate)) ** 2
    sel = (f >= 1.0) & (f <= 30.0)
    print("spectral peak (unconstrained 1-30 Hz): %.2f Hz" % f[sel][np.argmax(p[sel])])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-tc", type=int, default=10)
    ap.add_argument("--n-re", type=int, default=10)
    ap.add_argument("--tstop", type=float, default=12000.0)
    ap.add_argument("--gh-tc", type=float, default=4e-4)
    ap.add_argument("--gsk-re", type=float, default=5e-5)
    ap.add_argument("--out", default="")
    ap.add_argument("--bench", action="store_true",
                    help="scaling benchmark: 100/500/2000/5000 cells, 1 s each")
    a = ap.parse_args()

    if a.bench:
        pc = h.ParallelContext()
        rank, nhost = int(pc.id()), int(pc.nhost())
        if rank == 0:
            print("=== NEURON MPI scaling benchmark (%d ranks) ===" % nhost)
            print("replaces the laptop N^1.67 extrapolation with a real "
                  "measurement\n")
            print("%-9s %-9s %-11s %-13s" % ("cells", "wall (s)", "x realtime",
                                             "s per 1k cells"))
            print("-" * 46)
        for n in (100, 500, 2000, 5000):
            net = ParallelThalamicNet(n_tc=n // 2, n_re=n // 2,
                                      gh_tc=a.gh_tc, gsk_re=a.gsk_re)
            wall = net.run(tstop=1000.0)
            if rank == 0:
                print("%-9d %-9.1f %-11.2f %-13.2f"
                      % (n, wall, wall / 1.0, wall / (n / 1000.0)))
            del net
        if rank == 0:
            print("-" * 46)
        pc.barrier()
        pc.done()
        h.quit()          # required: without it MPI_Finalize is skipped -> abort
        return

    net = ParallelThalamicNet(n_tc=a.n_tc, n_re=a.n_re,
                              gh_tc=a.gh_tc, gsk_re=a.gsk_re)
    wall = net.run(tstop=a.tstop)
    t, g = net.gather()
    if net.rank == 0:
        _report(net, t, g, wall, a.tstop)
        if a.out:
            os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
            np.savez_compressed(a.out.replace(".h5", ".npz"),
                                times=t, gids=g, n_tc=net.n_tc, n_re=net.n_re,
                                tstop=a.tstop, wall=wall, gh_tc=a.gh_tc,
                                gsk_re=a.gsk_re)
            print("wrote %s" % a.out.replace(".h5", ".npz"))
    net.pc.barrier()
    net.pc.done()
    h.quit()          # required: without it MPI_Finalize is skipped -> abort


if __name__ == "__main__":
    main()
