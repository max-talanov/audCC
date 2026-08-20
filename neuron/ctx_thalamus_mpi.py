"""
MPI-parallel corticothalamic network for MareNostrum 5 -- the scale-out path
for ctx_thalamus_network.CorticoThalamicNet (10 TC + 10 RE + a ~40-cell
cortical column on one process) up to a bio-plausible column size.

Population sizes default to config/network_auditory_mn5_5k.yaml's DECLARED
sizes (3050 cells; the NEST reference model reaches ~5010 after its internal
RS/FRB/TuftRS-TuftIB/Basket-LTS-Axoaxonic subtype splits, which this module
does not reproduce 1:1 -- see --scale to go further):

    thalamus TC 210 / RE 55
    L4       E 640 / I 160
    L2/3     E 640 / I 160
    L5       E 530 / I 130   (half TuftIB via cortex_neuron.PYCellIB)
    L6       E 420 / I 105

Same pattern as tc_mpi.py (ParallelContext, gid-based registration, FIXED
convergence -- not the serial model's `frac`-based fan-in, which would give
every cell a fan-in of hundreds-to-thousands at this scale and average
discrete volleys into a tonic current; see tc_mpi.py's own note on this).
Cells are the SAME objects as the serial model (cortex_neuron.PYCell /
PYCellIB / FSCell, tc_neuron.TCCell / RECell), so the mechanism set is
identical: hh2, it/it2, ical, inap, cad, sk2, ihca.

Run:
    # scaling benchmark
    srun -n 100 .venv-neuron/bin/nrniv -python -mpi ctx_thalamus_mpi.py --bench

    # production at the reference-model size
    srun -n 100 .venv-neuron/bin/nrniv -python -mpi ctx_thalamus_mpi.py \
         --tstop 200000 --out out/mn5_ctx_nrn.h5

LIMITATION -- gap junctions (RE<->RE only) are rank-local, same caveat as
tc_mpi.py: see neuron/README.md.
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from neuron import h                                    # noqa: E402

h.nrnmpi_init()

import tc_neuron as T                                    # noqa: E402
import cortex_neuron as C                                 # noqa: E402

h.load_file("stdrun.hoc")

# Declared sizes, config/network_auditory_mn5_5k.yaml (thalamus + L4/L2-3/L5/L6).
DEFAULT_SIZES = {
    "tc": 210, "re": 55,
    "l4e": 640, "l4i": 160,
    "l23e": 640, "l23i": 160,
    "l5e": 530, "l5i": 130,
    "l6e": 420, "l6i": 105,
}


class ParallelCorticoThalamicNet:
    """Thalamic TC<->RE loop + a full cortical column (L4/L2-3/L5/L6),
    distributed over MPI ranks with fixed-convergence wiring.

    gid layout (contiguous blocks, in this order):
        tc, re, l4e, l4i, l23e, l23i, l5e, l5i, l6e, l6i
    """

    POPS = ["tc", "re", "l4e", "l4i", "l23e", "l23i", "l5e", "l5i", "l6e", "l6i"]

    def __init__(self, sizes=None, seed=1,
                 g_re_tc=0.015, g_tc_re=0.011, g_re_re=0.002, g_gap=0.03,
                 gh_tc=0.0, gsk_re=5e-5,
                 g_ff=0.0015, g_l5_l6=0.009, g_l5_l23=0.0015,
                 g_e_i=0.02, g_i_e=0.08, gsk_cx=8e-4, ib_frac=0.5,
                 g_tc_l4=0.0005, g_l6_tc=0.03, g_l6_re=0.03,
                 conv=100, gap_deg=6, gap_short=2,
                 het=0.05, delay_jitter=0.0):
        self.pc = h.ParallelContext()
        self.rank = int(self.pc.id())
        self.nhost = int(self.pc.nhost())
        self.sizes = dict(DEFAULT_SIZES)
        if sizes:
            self.sizes.update(sizes)
        self.rng = np.random.default_rng(seed + self.rank)
        self.cells, self.syns = {}, {}
        self.ncs, self.stims, self.gaps = [], [], []
        self.ib_frac = ib_frac
        # het/delay_jitter: SAME per-cell heterogeneity as the serial model
        # (ctx_thalamus_network.CorticoThalamicNet), which measured a 3x gain
        # in RE<->TC oscillatory cycles per SO event at het=0.05 (see
        # neuron/README.md "Heterogeneity"). Ported here so a production run
        # at bio-plausible scale reflects that finding rather than silently
        # regressing to the homogeneous (2.5 cycles/event) behaviour.
        # CRITICAL: jitter must be a DETERMINISTIC function of gid (its own
        # RNG stream, not self.rng), or the result depends on which rank
        # happens to build a given cell/connection -- exactly the rank-count
        # dependence the fixed-convergence design (see _draw's docstring
        # note) exists to avoid, and the property just verified on MN5
        # (identical spike counts at 2/4/100 ranks).
        self.het, self.delay_jitter = het, delay_jitter

        # -- gid ranges (contiguous blocks) ----------------------------------
        self.ranges = {}
        lo = 0
        for pop in self.POPS:
            n = self.sizes[pop]
            self.ranges[pop] = (lo, lo + n)
            lo += n
        self.n_total = lo

        # -- build cells: round-robin across ranks for load balance, EXCEPT
        #    re (block-assigned, so gap junctions cluster per rank, as in
        #    tc_mpi.py) --------------------------------------------------
        for pop in self.POPS:
            plo, phi = self.ranges[pop]
            n = phi - plo
            if pop == "re":
                per = int(np.ceil(n / self.nhost))
                mylo = plo + self.rank * per
                gids = list(range(mylo, min(mylo + per, phi)))
            else:
                gids = list(range(plo + self.rank, phi, self.nhost))
            for gid in gids:
                c = self._make_cell(pop, gid)
                self._register(gid, c)

        # -- synapses + fixed-convergence wiring -----------------------------
        self.conv = conv
        self._wire_re_tc(g_re_tc)
        self._wire_tc_re(g_tc_re)
        self._wire_re_re_local(g_re_re)
        self._wire_gap(g_gap, gap_deg, gap_short)
        self._wire_intracortical(g_ff, g_l5_l6, g_l5_l23)
        for pop_e, pop_i in [("l4e", "l4i"), ("l23e", "l23i"),
                              ("l5e", "l5i"), ("l6e", "l6i")]:
            self._wire_layer_inh(pop_e, pop_i, g_e_i, g_i_e)
        self._wire_thalamocortical(g_tc_l4)
        self._wire_corticothalamic(g_l6_tc, g_l6_re)

        self.tspk, self.gspk = h.Vector(), h.Vector()
        self.pc.spike_record(-1, self.tspk, self.gspk)

    # ------------------------------------------------------------- cell build
    JITTER_OFFSET = 10_000_000  # keeps _jitter's per-gid RNG stream disjoint
                                 # from _draw's per-gid connectivity streams

    def _jitter(self, base, gid, salt=0):
        """Deterministic function of (gid, salt) -- NOT self.rng -- so the
        result is independent of which rank built this cell/connection."""
        if self.het == 0:
            return base
        r = np.random.default_rng(self.JITTER_OFFSET + gid * 97 + salt)
        return base * (1.0 + r.uniform(-self.het, self.het))

    def _make_cell(self, pop, gid):
        if pop == "tc":
            c = T.TCCell(gsk=0.0, gh=0.0)
            c.soma.e_pas = self._jitter(-80.0, gid, 1)
            c.soma.gcabar_it *= self._jitter(1.0, gid, 2)
            return c
        if pop == "re":
            c = T.RECell(gsk=5e-5)
            c.soma.e_pas = self._jitter(-82.0, gid, 1)
            c.soma.gcabar_it2 *= self._jitter(1.0, gid, 2)
            return c
        if pop.endswith("i"):
            return C.FSCell(e_pas=self._jitter(-68.0, gid, 1))
        # excitatory cortical: L5 gets a deterministic IB fraction (first
        # ib_frac of each layer's local gid range) so it is reproducible
        # across rank counts, matching the fixed-convergence philosophy.
        plo, phi = self.ranges[pop]
        n = phi - plo
        n_ib = int(round(n * self.ib_frac)) if pop == "l5e" else 0
        if gid - plo < n_ib:
            return C.PYCellIB(e_pas=self._jitter(-70.0, gid, 1),
                               gnap=self._jitter(2e-4, gid, 3))
        return C.PYCell(e_pas=self._jitter(-70.0, gid, 1), gsk=8e-4)

    # ------------------------------------------------------------------ utils
    def _register(self, gid, cell):
        self.cells[gid] = cell
        self.pc.set_gid2node(gid, self.rank)
        nc = h.NetCon(cell.soma(0.5)._ref_v, None, sec=cell.soma)
        nc.threshold = -10
        self.pc.cell(gid, nc)
        self.ncs.append(nc)

    def _target_sec(self, gid):
        tgt = getattr(self.cells[gid], "synsec", None)
        return tgt if tgt is not None else self.cells[gid].soma

    def _draw(self, seed_gid, lo, hi, k):
        r = np.random.default_rng(seed_gid)
        n = hi - lo
        k = min(k, n)
        return (lo + r.choice(n, size=k, replace=False)) if n and k else []

    def _connect(self, src_gid, syn, weight, delay, dst_gid=None):
        if self.delay_jitter > 0 and dst_gid is not None:
            r = np.random.default_rng(self.JITTER_OFFSET * 2
                                       + int(src_gid) * 97 + int(dst_gid))
            delay = max(0.1, delay + r.uniform(-self.delay_jitter,
                                                self.delay_jitter))
        nc = self.pc.gid_connect(int(src_gid), syn)
        nc.weight[0], nc.delay = weight, delay
        self.ncs.append(nc)

    def _project(self, key, pre_pop, post_pop, e, tau1, tau2, g, delay=1.0,
                 seed_offset=0):
        """Fixed-convergence excitatory/inhibitory projection, cells this
        rank owns as TARGET only (NEURON connects TO a local gid)."""
        plo, phi = self.ranges[pre_pop]
        qlo, qhi = self.ranges[post_pop]
        k = max(1, min(self.conv, phi - plo))
        for gid in range(qlo, qhi):
            if int(self.pc.gid_exists(gid)) == 0:
                continue
            syn = h.Exp2Syn(self._target_sec(gid)(0.5))
            syn.e, syn.tau1, syn.tau2 = e, tau1, tau2
            self.syns[(key, gid)] = syn
            for src in self._draw(gid + seed_offset, plo, phi, k):
                self._connect(src, syn, g / k, delay, dst_gid=gid)

    # -- intrathalamic wiring (params matched to ctx_thalamus_network.py) --
    def _wire_re_tc(self, g):
        self._project("gabaa", "re", "tc", -85.0, 1.0, 8.0, g, seed_offset=0)

    def _wire_tc_re(self, g):
        self._project("ampa", "tc", "re", 0.0, 0.5, 2.0, g, seed_offset=1)

    def _wire_re_re_local(self, g):
        """Nearest-neighbour lateral inhibition (|i-j|<=2), fixed degree --
        NOT convergence-scaled, matching tc_network_nrn._wire_re_re."""
        rlo, rhi = self.ranges["re"]
        n_re = rhi - rlo
        for gid in range(rlo, rhi):
            if int(self.pc.gid_exists(gid)) == 0:
                continue
            syn = h.Exp2Syn(self.cells[gid].soma(0.5))
            syn.e, syn.tau1, syn.tau2 = -75.0, 1.0, 6.0
            self.syns[("gabaa_re", gid)] = syn
            i = gid - rlo
            for j in range(max(0, i - 2), min(n_re, i + 3)):
                if j != i:
                    self._connect(rlo + j, syn, g, 1.0, dst_gid=gid)

    def _wire_gap(self, g, gap_deg, gap_short):
        """Cross-rank gap junctions between RE cells (pc.setup_transfer),
        small-world topology -- see tc_mpi.py for the full rationale."""
        rlo, rhi = self.ranges["re"]
        n_re = rhi - rlo
        if g <= 0 or n_re <= 1:
            return
        SGID = self.n_total + 10000
        for gid in range(rlo, rhi):
            if int(self.pc.gid_exists(gid)) == 0:
                continue
            sec = self.cells[gid].soma
            self.pc.source_var(sec(0.5)._ref_v, SGID + (gid - rlo), sec=sec)
        half = max(1, gap_deg // 2)
        for gid in range(rlo, rhi):
            if int(self.pc.gid_exists(gid)) == 0:
                continue
            i = gid - rlo
            nbrs = set()
            for d in range(1, half + 1):
                nbrs.add((i + d) % n_re)
                nbrs.add((i - d) % n_re)
            if gap_short > 0 and n_re > 2 * half + 1:
                r = np.random.default_rng(2000003 + i)
                for _ in range(gap_short):
                    nbrs.add(int(r.integers(0, n_re)))
            nbrs.discard(i)
            for nb in nbrs:
                gap = h.GapMPI(self.cells[gid].soma(0.5))
                gap.g = g
                self.pc.target_var(gap, gap._ref_vgap, SGID + nb)
                self.gaps.append(gap)
        self.pc.setup_transfer()

    # -- intracortical wiring (params matched to cortex_neuron.CorticalColumn) --
    def _wire_intracortical(self, g_ff, g_l5_l6, g_l5_l23):
        self._project("l4_l23", "l4e", "l23e", 0.0, 0.5, 2.0, g_ff, seed_offset=10)
        self._project("l23_l5", "l23e", "l5e", 0.0, 0.5, 2.0, g_ff, seed_offset=11)
        self._project("l4_l5", "l4e", "l5e", 0.0, 0.5, 2.0, g_ff * 0.5, seed_offset=12)
        self._project("l5_l6", "l5e", "l6e", 0.0, 0.5, 2.0, g_l5_l6, seed_offset=13)
        self._project("l5_l23", "l5e", "l23e", 0.0, 0.5, 2.0, g_l5_l23, seed_offset=14)

    def _wire_layer_inh(self, pop_e, pop_i, g_e_i, g_i_e):
        self._project(f"{pop_e}_i", pop_e, pop_i, 0.0, 0.5, 2.0, g_e_i, seed_offset=20)
        self._project(f"{pop_i}_e", pop_i, pop_e, -75.0, 1.0, 8.0, g_i_e, seed_offset=21)

    # -- thalamocortical / corticothalamic loop --------------------------
    def _wire_thalamocortical(self, g):
        self._project("tc_l4", "tc", "l4e", 0.0, 0.5, 2.0, g, seed_offset=30)

    def _wire_corticothalamic(self, g_tc, g_re):
        self._project("l6_tc", "l6e", "tc", 0.0, 0.5, 2.0, g_tc, seed_offset=31)
        self._project("l6_re", "l6e", "re", 0.0, 0.5, 2.0, g_re, seed_offset=32)

    # -------------------------------------------------------------------- run
    def run(self, tstop=12000.0, celsius=36.0, dt=0.025):
        self.pc.set_maxstep(10)
        h.celsius, h.dt = celsius, dt
        h.finitialize(-70)
        t0 = time.time()
        self.pc.psolve(tstop)
        wall = time.time() - t0
        wall = self.pc.allreduce(wall, 2)
        return wall

    def teardown(self):
        self.cells.clear()
        self.syns.clear()
        self.ncs.clear()
        self.stims.clear()
        self.gaps.clear()
        self.tspk = self.gspk = None
        import gc
        gc.collect()
        self.pc.gid_clear()

    def gather(self):
        t_all = self.pc.py_gather(np.asarray(self.tspk), 0)
        g_all = self.pc.py_gather(np.asarray(self.gspk), 0)
        if self.rank != 0:
            return None, None
        t = np.concatenate(t_all) if t_all else np.array([])
        g = np.concatenate(g_all) if g_all else np.array([])
        o = np.argsort(t)
        return t[o], g[o]


def _report(net, t, g, wall, tstop):
    print("\n=== NEURON MPI corticothalamic network ===")
    print("ranks %d | total cells %d | tstop %.0f ms | wall %.1f s (%.2f x realtime)"
          % (net.nhost, net.n_total, tstop, wall, wall / (tstop / 1000.0)))
    for pop in net.POPS:
        lo, hi = net.ranges[pop]
        print("  %-6s gid [%d,%d)  n=%d" % (pop, lo, hi, hi - lo))
    if t is None or len(t) < 20:
        print("too few spikes to analyse (%d)" % (0 if t is None else len(t)))
        return
    m = t > min(2000.0, 0.2 * tstop)  # discard startup transient (relative,
    # so a short smoke-test tstop doesn't get its whole spike train discarded)
    t, g = t[m], g[m]
    tc_lo, tc_hi = net.ranges["tc"]
    re_lo, re_hi = net.ranges["re"]
    l6_lo, l6_hi = net.ranges["l6e"]
    n_tc_sp = int(((g >= tc_lo) & (g < tc_hi)).sum())
    n_re_sp = int(((g >= re_lo) & (g < re_hi)).sum())
    n_l6_sp = int(((g >= l6_lo) & (g < l6_hi)).sum())
    print("spikes %d total (TC %d, RE %d, L6E %d)" % (len(t), n_tc_sp, n_re_sp, n_l6_sp))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tstop", type=float, default=12000.0)
    ap.add_argument("--conv", type=int, default=100)
    ap.add_argument("--scale", type=float, default=1.0,
                    help="multiply every population size (default sizes sum "
                         "to the reference model's declared 3050)")
    ap.add_argument("--out", default="")
    ap.add_argument("--bench", action="store_true",
                    help="scaling benchmark: several sizes, 1 s each")
    ap.add_argument("--bench-scales", default="0.1,0.5,1.0,1.65",
                    help="comma-separated --scale values for --bench "
                         "(1.65x default ~= 5010, the NEST model's actual size)")
    ap.add_argument("--bench-ms", type=float, default=1000.0)
    ap.add_argument("--het", type=float, default=0.05,
                    help="per-cell heterogeneity, uniform +/- fraction on "
                         "e_pas/gcabar_it(2)/gnap (default matches the "
                         "serial model's validated 0.05; see neuron/README.md "
                         "'Heterogeneity'). Deterministic per-gid, so results "
                         "stay rank-count independent.")
    ap.add_argument("--delay-jitter", type=float, default=0.0,
                    help="+/- ms uniform jitter on every synaptic delay")
    a = ap.parse_args()

    if a.bench:
        pc = h.ParallelContext()
        rank, nhost = int(pc.id()), int(pc.nhost())
        if rank == 0:
            print("=== NEURON MPI corticothalamic scaling benchmark (%d ranks) ==="
                  % nhost)
            # spikes/gid breakdown is printed alongside timing SPECIFICALLY so
            # this doubles as the rank-count correctness check: run the same
            # --bench-scales at --ntasks 2 and --ntasks 4 and diff the spike
            # counts, not just the timing (timing alone says nothing about
            # whether cross-rank gid_connect / gap junctions are wired right).
            print("%-9s %-9s %-11s %-13s %-9s %-6s %-6s %-6s"
                  % ("cells", "wall (s)", "x realtime", "s per 1k cells",
                     "spikes", "TC", "RE", "L6E"))
            print("-" * 76)
        for s in [float(x) for x in a.bench_scales.split(",")]:
            sizes = {k: max(1, int(round(v * s))) for k, v in DEFAULT_SIZES.items()}
            net = ParallelCorticoThalamicNet(sizes=sizes, conv=a.conv,
                                              het=a.het,
                                              delay_jitter=a.delay_jitter)
            wall = net.run(tstop=a.bench_ms)
            t, g = net.gather()
            if rank == 0:
                tc_lo, tc_hi = net.ranges["tc"]
                re_lo, re_hi = net.ranges["re"]
                l6_lo, l6_hi = net.ranges["l6e"]
                n = len(t) if t is not None else 0
                n_tc = int(((g >= tc_lo) & (g < tc_hi)).sum()) if n else 0
                n_re = int(((g >= re_lo) & (g < re_hi)).sum()) if n else 0
                n_l6 = int(((g >= l6_lo) & (g < l6_hi)).sum()) if n else 0
                print("%-9d %-9.1f %-11.2f %-13.2f %-9d %-6d %-6d %-6d"
                      % (net.n_total, wall, wall / (a.bench_ms / 1000.0),
                         wall / (net.n_total / 1000.0), n, n_tc, n_re, n_l6))
            net.teardown()
            del net
        if rank == 0:
            print("-" * 76)
        pc.barrier()
        pc.done()
        h.quit()
        return

    sizes = {k: max(1, int(round(v * a.scale))) for k, v in DEFAULT_SIZES.items()}
    net = ParallelCorticoThalamicNet(sizes=sizes, conv=a.conv, het=a.het,
                                      delay_jitter=a.delay_jitter)
    wall = net.run(tstop=a.tstop)
    t, g = net.gather()
    if net.rank == 0:
        _report(net, t, g, wall, a.tstop)
        if a.out:
            os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
            np.savez_compressed(a.out.replace(".h5", ".npz"),
                                 times=t, gids=g, sizes=net.sizes,
                                 ranges=net.ranges, tstop=a.tstop, wall=wall)
    net.teardown()
    pc = h.ParallelContext()
    pc.barrier()
    pc.done()
    h.quit()


if __name__ == "__main__":
    main()
