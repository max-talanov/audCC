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
                 g_re_tc=0.015, tau2_re_tc=8.0, g_tc_re=0.011, g_re_re=0.006, g_re_re_sd=0.002,
                 g_gap=0.03,
                 gh_tc=0.0, gsk_re=1e-3,
                 g_ff=0.0015, g_l5_l6=0.009, g_l5_l23=0.0015,
                 g_l5_rec=0.0, tau2_l5_rec=120.0, l5_rec_mech="exp2syn",
                 mg_l5_rec=1.0, tau1_l5_rec_nmda=5.0,
                 g_e_i=0.02, g_i_e=0.08, g_i_e_l5=0.0, gsk_cx=8e-4, ib_frac=0.5,
                 g_tc_l4=0.02, g_l6_tc=0.03, g_l6_re=0.03,
                 conv=100, gap_deg=6, gap_short=2, g_l5_gap=0.02,
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
        # gsk_re was DECLARED but never used -- _make_cell hardcoded 5e-5, a
        # value fitted against the old GHK itd port. On the published (ohmic)
        # it2 that gives 20-spike bursts @ 317 Hz, and with L6 feedback driving
        # them the reticular population ran away to 230 Hz/cell for the whole
        # 200 s of job 44895771 (78% of all spikes from 91 cells). Refitted
        # value is 1e-3 -> 8 spikes @ 363 Hz.
        self.gsk_re = gsk_re
        # gh_tc was ALSO declared but never used (same bug as gsk_re above):
        # _make_cell hardcoded gh=0.0, so TC's Ca2+-dependent I_h (ihca.mod)
        # has never actually been enabled in the MPI/production model despite
        # the CLI/constructor accepting a value for it.
        self.gh_tc = gh_tc
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
        self._gap_wired = False
        self._wire_re_tc(g_re_tc, tau2_re_tc)
        self._wire_tc_re(g_tc_re)
        self._wire_re_re_local(g_re_re, g_re_re_sd)
        self._wire_gap(g_gap, gap_deg, gap_short)
        self._wire_l5_gap(g_l5_gap, gap_deg, gap_short)
        if self._gap_wired:
            # ONE setup_transfer() call after ALL source_var/target_var
            # registrations (RE and L5 gaps both), not one per group.
            self.pc.setup_transfer()
        self._wire_intracortical(g_ff, g_l5_l6, g_l5_l23)
        self._wire_l5_recurrent(g_l5_rec, tau2_l5_rec, mech=l5_rec_mech,
                                 mg=mg_l5_rec, tau1_nmda=tau1_l5_rec_nmda)
        for pop_e, pop_i in [("l4e", "l4i"), ("l23e", "l23i"),
                              ("l5e", "l5i"), ("l6e", "l6i")]:
            # L5's own I->E feedback (g_i_e_l5) is kept separate from the
            # other layers' (g_i_e): a small-scale isolation test found that
            # the shared g_i_e=0.08 recruits L5I (FSCell) repeatedly off the
            # near-synchronous IB population's first 1-2 spikes, holding
            # L5E clamped near the GABA_A reversal (-75 mV) for tens of ms --
            # long enough to defeat I_NaP before the plateau burst (~5-8
            # spikes, ~10-500 ms, the article's own slow-oscillation UP
            # state) can develop at all, and before the SLOWED SK2/Ca2+
            # pool (taur=500 ms) gets a chance to be what actually
            # terminates it, as PYCellIB's docstring intends. With L5's
            # feedback inhibition removed the plateau is restored (5-6
            # spikes, ~11-13 ms span, clean ~0.77 Hz UP/DOWN cycle, stable
            # over 6+ cycles/8s) -- matching the isolated single-cell
            # behaviour almost exactly. g_i_e_l5 stays a separate knob
            # (not simply "g_i_e=0 everywhere") since L4/L2-3/L6 don't have
            # this intrinsic-bursting mechanism and still need their own
            # feedback inhibition for stability.
            gie = g_i_e_l5 if pop_e == "l5e" else g_i_e
            self._wire_layer_inh(pop_e, pop_i, g_e_i, gie)
        self._wire_thalamocortical(g_tc_l4)
        self._wire_corticothalamic(g_l6_tc, g_l6_re)

        self.tspk, self.gspk = h.Vector(), h.Vector()
        self.pc.spike_record(-1, self.tspk, self.gspk)

    # ------------------------------------------------------------- cell build
    JITTER_OFFSET = 10_000_000  # keeps _jitter's per-gid RNG stream disjoint
                                 # from _draw's per-gid connectivity streams
    WEIGHT_JITTER_OFFSET = 30_000_000  # per-connection synaptic weight jitter
                                         # (distinct from JITTER_OFFSET*2, the
                                         # delay-jitter stream in _connect)

    def _jitter(self, base, gid, salt=0):
        """Deterministic function of (gid, salt) -- NOT self.rng -- so the
        result is independent of which rank built this cell/connection."""
        if self.het == 0:
            return base
        r = np.random.default_rng(self.JITTER_OFFSET + gid * 97 + salt)
        return base * (1.0 + r.uniform(-self.het, self.het))

    def _make_cell(self, pop, gid):
        if pop == "tc":
            c = T.TCCell(gsk=0.0, gh=self.gh_tc)
            c.soma.e_pas = self._jitter(-80.0, gid, 1)
            c.soma.gcabar_it *= self._jitter(1.0, gid, 2)
            return c
        if pop == "re":
            c = T.RECell(gsk=self.gsk_re)
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
            # PYCellIB gets the same jitter as every other cell type --
            # individual variability is the biologically correct default
            # (real neurons are not copies of one ideal cell). An A/B test
            # (scale=0.3, 159 L5E cells) found that het on e_pas/gnap alone
            # desynchronises the IB population once it's a few hundred
            # cells (L6E/RE went to ~0 spikes on the full 5031-cell MN5 run,
            # job 44844450) -- but the fix for that is an explicit
            # synchronising mechanism (gap junctions between IB cells, see
            # _wire_l5_gap below), the same way TRN heterogeneity + RE<->RE
            # gap junctions coexist in tc_neuron.py, not making the cells
            # identical.
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
    def _wire_re_tc(self, g, tau2=8.0):
        # tau2: RE->TC GABA_A decay. The IPSP that hyperpolarises TC long
        # enough to de-inactivate I_T for a rebound burst -- candidate lever
        # for RE/TC multi-cycle ringing (a longer IPSP -> a longer, possibly
        # more effective de-inactivation window before rebound) now that
        # gh_tc was ruled out as that lever (see --sweep-gh-tc).
        self._project("gabaa", "re", "tc", -85.0, 1.0, tau2, g, seed_offset=0)

    def _wire_tc_re(self, g):
        self._project("ampa", "tc", "re", 0.0, 0.5, 2.0, g, seed_offset=1)

    def _wire_re_re_local(self, g, g_sd=0.0):
        """Nearest-neighbour lateral inhibition (|i-j|<=2), fixed degree --
        NOT convergence-scaled, matching tc_network_nrn._wire_re_re.

        g_sd > 0: each individual RE->RE synapse's weight is drawn from
        N(g, g_sd) instead of every connection sharing the identical value g.
        Real synapses are not identical -- and unlike jittering an intrinsic
        cell parameter (the L5 IB pacemaker lesson), this is a MEAN-preserving
        per-connection distribution: each RE cell still sums inhibition from
        ~4 neighbours, so the aggregate a cell receives stays close to the
        mean even though individual connections vary. Deterministic per
        (src, dst) pair via its own RNG stream (WEIGHT_JITTER_OFFSET), so the
        result is independent of which rank builds the connection -- same
        requirement as _jitter/_connect's delay jitter. Truncated at a small
        positive floor so no connection goes to zero or negative."""
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
                    src = rlo + j
                    w = g
                    if g_sd > 0:
                        r = np.random.default_rng(self.WEIGHT_JITTER_OFFSET
                                                   + int(src) * 97 + int(gid))
                        w = max(g * 0.05, r.normal(g, g_sd))
                    self._connect(src, syn, w, 1.0, dst_gid=gid)

    def _wire_gap(self, g, gap_deg, gap_short):
        """Cross-rank gap junctions between RE cells (pc.setup_transfer),
        small-world topology -- see tc_mpi.py for the full rationale."""
        rlo, rhi = self.ranges["re"]
        self._wire_gap_range(rlo, rhi, g, gap_deg, gap_short,
                              sgid_base=self.n_total + 10000, seed_base=2000003)

    def _wire_l5_gap(self, g, gap_deg, gap_short):
        """Cross-rank gap junctions between L5's intrinsically-bursting
        cells -- the explicit synchronising mechanism that lets the IB
        population stay phase-locked despite per-cell heterogeneity in its
        own burst-timing parameters (e_pas, gnap). Same small-world topology
        as _wire_gap; without this, the IB population desynchronises once
        it's a few hundred cells and the L5->L6->TC/RE loop goes dead (see
        neuron/README.md 'MN5 5k-scale result')."""
        lo, hi = self.ranges["l5e"]
        n_ib = int(round((hi - lo) * self.ib_frac))
        self._wire_gap_range(lo, lo + n_ib, g, gap_deg, gap_short,
                              sgid_base=self.n_total + 20000, seed_base=3000007)

    def _wire_gap_range(self, lo, hi, g, gap_deg, gap_short, sgid_base, seed_base):
        n = hi - lo
        if g <= 0 or n <= 1:
            return
        for gid in range(lo, hi):
            if int(self.pc.gid_exists(gid)) == 0:
                continue
            sec = self.cells[gid].soma
            self.pc.source_var(sec(0.5)._ref_v, sgid_base + (gid - lo), sec=sec)
        half = max(1, gap_deg // 2)
        for gid in range(lo, hi):
            if int(self.pc.gid_exists(gid)) == 0:
                continue
            i = gid - lo
            nbrs = set()
            for d in range(1, half + 1):
                nbrs.add((i + d) % n)
                nbrs.add((i - d) % n)
            if gap_short > 0 and n > 2 * half + 1:
                r = np.random.default_rng(seed_base + i)
                for _ in range(gap_short):
                    nbrs.add(int(r.integers(0, n)))
            nbrs.discard(i)
            for nb in nbrs:
                gap = h.GapMPI(self.cells[gid].soma(0.5))
                gap.g = g
                self.pc.target_var(gap, gap._ref_vgap, sgid_base + nb)
                self.gaps.append(gap)
        self._gap_wired = True

    # -- intracortical wiring (params matched to cortex_neuron.CorticalColumn) --
    def _wire_intracortical(self, g_ff, g_l5_l6, g_l5_l23):
        self._project("l4_l23", "l4e", "l23e", 0.0, 0.5, 2.0, g_ff, seed_offset=10)
        self._project("l23_l5", "l23e", "l5e", 0.0, 0.5, 2.0, g_ff, seed_offset=11)
        self._project("l4_l5", "l4e", "l5e", 0.0, 0.5, 2.0, g_ff * 0.5, seed_offset=12)
        self._project("l5_l6", "l5e", "l6e", 0.0, 0.5, 2.0, g_l5_l6, seed_offset=13)
        self._project("l5_l23", "l5e", "l23e", 0.0, 0.5, 2.0, g_l5_l23, seed_offset=14)

    REC_OFFSET = 40_000_000  # per-gid RNG stream for recurrent wiring,
                              # disjoint from JITTER_OFFSET/WEIGHT_JITTER_OFFSET

    def _wire_l5_recurrent(self, g, tau2=120.0, conv=None, mech="exp2syn",
                            mg=1.0, tau1_nmda=5.0):
        """Recurrent L5E->L5E excitation -- test of the "missing slow
        recurrent excitation" hypothesis for a genuine cortical UP state
        (see architecture review: no layer has ANY within-population E->E
        synapse, and no synapse anywhere is slower than 8ms -- so nothing
        can sustain elevated firing once triggered; L5's I_NaP plateau is a
        lone single-cell/gap-junction pacemaker with no network-level
        reinforcement, which is why a single fast IPSP was enough to kill
        it outright before the g_i_e_l5 fix).

        mech="exp2syn" (Option B): approximates NMDA-like sustaining
        excitation with the SAME Exp2Syn used everywhere else, just with a
        much longer tau2 (default 120ms vs. 2-8ms elsewhere) -- cheap,
        no new mechanism. Result: a real, powerful lever, but with
        essentially NO stable middle ground -- g_l5_rec<=0.0005 barely
        changes anything and >=0.001 causes runaway/continuous tonic
        firing at every g_i_e_l5 tested (0.01-0.08). A plain linear
        conductance has no self-limiting nonlinearity, so nothing caps the
        positive feedback except inhibition racing to catch up.

        mech="nmda" (Option A): a real NMDA mechanism (mod/nmda.mod,
        Jahr & Stevens 1990 Mg2+ block) -- voltage-dependent conductance,
        weak near rest and unblocking with depolarisation, which is
        exactly the saturating nonlinearity a STABLE bistable UP state
        needs and Option B's result argues is functionally necessary, not
        just a bio-plausibility upgrade.

        g=0.0 (the default) leaves this fully disabled -- opt-in only.
        """
        if g <= 0:
            return
        lo, hi = self.ranges["l5e"]
        n = hi - lo
        k = max(1, min(conv or self.conv, n - 1))
        for gid in range(lo, hi):
            if int(self.pc.gid_exists(gid)) == 0:
                continue
            if mech == "nmda":
                syn = h.NMDA(self._target_sec(gid)(0.5))
                syn.tau1, syn.tau2, syn.mg = tau1_nmda, tau2, mg
            else:
                syn = h.Exp2Syn(self._target_sec(gid)(0.5))
                syn.e, syn.tau1, syn.tau2 = 0.0, 0.5, tau2
            self.syns[("l5_rec", gid)] = syn
            r = np.random.default_rng(self.REC_OFFSET + gid)
            others = np.array([x for x in range(lo, hi) if x != gid])
            srcs = r.choice(others, size=min(k, len(others)), replace=False)
            for src in srcs:
                self._connect(int(src), syn, g / k, 1.0, dst_gid=gid)

    def _wire_layer_inh(self, pop_e, pop_i, g_e_i, g_i_e):
        self._project(f"{pop_e}_i", pop_e, pop_i, 0.0, 0.5, 2.0, g_e_i, seed_offset=20)
        self._project(f"{pop_i}_e", pop_i, pop_e, -75.0, 1.0, 8.0, g_i_e, seed_offset=21)

    # -- thalamocortical / corticothalamic loop --------------------------
    def _wire_thalamocortical(self, g):
        """TC -> L4E, the ASCENDING arm of the corticothalamic loop.

        g was 0.0005 -- the smallest conductance in the model, 60x below
        g_l6_tc and 3x below g_ff, then divided by conv. Since tc_l4 is L4E's
        ONLY excitatory input (L4E is a source in l4_l23/l4_l5 but never a
        target elsewhere), L4 and L2/3 were completely SILENT on job 44895771:
        2640 of 5031 cells produced zero spikes. The cortex ran purely off
        L5's intrinsically-bursting cells via l5_l6, so the loop that executed
        was L5(IB) -> L6 -> thalamus, with the thalamocortical arm dead -- a
        descending oscillator, not a closed loop. Measured at scale 0.12:
        L4E 0.00 Hz at 0.0005, 4.77 Hz at 0.02."""
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


def _re_burst_stats(t, g, re_lo, re_hi, burst_gap=30.0, event_gap=300.0):
    """RE burst-SHAPE diagnostics -- the article's SK2 mechanism terminates a
    multi-spike burst; a network that fires one isolated spike per cycle
    never exercises it at all. Reports the fraction of spikes that are part
    of a burst (ISI < burst_gap to the previous spike from the SAME cell),
    the mean burst size among cells that ever burst, and the population
    event structure (SO-cycle volleys)."""
    m = (g >= re_lo) & (g < re_hi)
    tp, gp = t[m], g[m]
    if len(tp) < 2:
        return dict(n_spikes=len(tp), frac_burst=0.0, mean_burst_size=0.0,
                    n_events=0, event_hz=0.0)
    order = np.argsort(gp, kind="stable")
    tp, gp = tp[order], gp[order]
    same_cell = np.diff(gp) == 0
    isi = np.diff(tp)
    in_burst = same_cell & (isi < burst_gap)
    frac_burst = float(in_burst.sum()) / len(tp)
    # mean burst size: count runs of consecutive in_burst=True per cell,
    # +1 for the spike that starts each run
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
    return dict(n_spikes=len(tp), frac_burst=round(frac_burst, 3),
                mean_burst_size=round(mean_burst, 2), n_events=len(starts),
                event_hz=round(event_hz, 3))


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
    ap.add_argument("--sweep-g-re-re", default="",
                    help="comma-separated g_re_re values to sweep, AT FULL "
                         "PRODUCTION SCALE (--scale), one run each, reporting "
                         "RE burst-shape stats (frac_burst, mean_burst_size, "
                         "event_hz) -- not just total spike count. Use a short "
                         "--sweep-tstop, not the full 200 s production length.")
    ap.add_argument("--sweep-tstop", type=float, default=20000.0,
                    help="tstop per sweep point (default 20 s -- enough SO "
                         "cycles for burst statistics without full-length cost)")
    ap.add_argument("--g-re-re", type=float, default=0.006,
                    help="mean RE<->RE lateral GABA_A weight. Sweep found "
                         "0.005-0.03 gives real 3-7 spike RE bursts (article's "
                         "'2 to >10'); 0.002 causes runaway tonic firing (229 "
                         "Hz), 0.1 over-suppresses to single spikes (0 bursts). "
                         "Default 0.006 sits between the two best-tested points "
                         "(0.005, 0.01: 6.86/6.20 spikes/burst).")
    ap.add_argument("--g-re-re-sd", type=float, default=0.002,
                    help="per-connection normal-distribution SD around "
                         "--g-re-re (0 = every RE<->RE synapse identical). "
                         "Default +/-0.002 makes 0.002-0.01 approx +/-2 SD -- "
                         "real synapses are not identical, and only individual "
                         "connections reach the low tail while each RE cell "
                         "still sums input from ~4 neighbours.")
    ap.add_argument("--g-i-e-l5", type=float, default=0.0,
                    help="L5's own I->E feedback weight (FSCell->PYCellIB/E), "
                         "kept separate from --g-i-e everywhere else. A "
                         "small-scale isolation test found the shared "
                         "g_i_e=0.08 recruits L5I repeatedly off the "
                         "near-synchronous IB population's first 1-2 spikes, "
                         "clamping L5E near GABA_A reversal for tens of ms and "
                         "defeating the I_NaP plateau (the article's slow- "
                         "oscillation UP-state generator) before it can "
                         "develop. Default 0.0 lets SK2/Ca2+ (taur=500 ms) be "
                         "what terminates the plateau, as PYCellIB intends --"
                         " restores a clean 5-6 spike, ~11-13 ms burst at a "
                         "stable ~0.77 Hz, matching the isolated single-cell "
                         "behaviour.")
    ap.add_argument("--sweep-gh-tc", default="",
                    help="comma-separated gh_tc values to sweep (full network, "
                         "--scale), reporting BOTH TC and RE burst-shape stats. "
                         "The values validated on the old thalamus-only "
                         "ThalamicNet (2e-5, 4e-4) both broke in this network "
                         "(TC tonic runaway / outright explosion) -- this "
                         "sweep is for finding a much smaller working point, "
                         "e.g. \"0,1e-6,3e-6,1e-5\". Use a short --sweep-tstop: "
                         "gh_tc>0 was measured up to 17x slower to simulate "
                         "than gh_tc=0 (975s wall for 20s simulated at "
                         "gh_tc=2e-5, --scale 0.05), likely from the ~3x spike "
                         "count feeding more NetCon delivery events, not the "
                         "integrator step itself.")
    ap.add_argument("--sweep-g-l5-rec", default="",
                    help="comma-separated L5E->L5E recurrent-excitation "
                         "weights to sweep (full network, --scale), reporting "
                         "L5E burst-shape stats (does the plateau actually "
                         "SUSTAIN now?) plus RE/TC for side effects. Option B "
                         "test of the missing-recurrent-excitation hypothesis "
                         "(see architecture review) -- approximates NMDA-like "
                         "sustaining excitation with existing Exp2Syn at a "
                         "much longer tau2 (--tau2-l5-rec) than any other "
                         "synapse in the model. 0 = fully disabled (default). "
                         "Try e.g. \"0,0.0005,0.001,0.002,0.005\".")
    ap.add_argument("--g-l5-rec", type=float, default=0.0,
                    help="L5E->L5E recurrent excitation weight. See "
                         "--sweep-g-l5-rec. 0.0 (default) disables it.")
    ap.add_argument("--tau2-l5-rec", type=float, default=120.0,
                    help="L5E->L5E recurrent excitation decay (ms) -- the "
                         "NMDA-like slow time constant this synapse needs to "
                         "sustain a real UP state (every other synapse in "
                         "the model decays in 2-8ms).")
    ap.add_argument("--l5-rec-mech", choices=["exp2syn", "nmda"], default="exp2syn",
                    help="L5E->L5E recurrent synapse mechanism. 'exp2syn' "
                         "(Option B, default): plain long-tau2 Exp2Syn, found "
                         "to have essentially no stable middle ground (no "
                         "effect or runaway). 'nmda' (Option A): real "
                         "Mg2+-block voltage-dependent NMDA (mod/nmda.mod) "
                         "-- the saturating nonlinearity a stable bistable "
                         "UP state needs.")
    ap.add_argument("--mg-l5-rec", type=float, default=1.0,
                    help="extracellular Mg2+ (mM) for --l5-rec-mech=nmda "
                         "(1 mM is physiological; higher = stronger block).")
    ap.add_argument("--sweep-tau2-re-tc", default="",
                    help="comma-separated RE->TC GABA_A decay (tau2, ms) "
                         "values to sweep (full network, --scale), reporting "
                         "TC and RE burst-shape stats. gh_tc (I_h) was ruled "
                         "out as the multi-cycle-ringing lever (--sweep-gh-tc: "
                         "monotonically bad, never helps); this tests whether "
                         "a longer RE->TC IPSP (more time to de-inactivate "
                         "I_T before rebound) does better. Default 8.0 ms; "
                         "try e.g. \"8,15,25,40,60,100\".")
    ap.add_argument("--tau2-re-tc", type=float, default=8.0,
                    help="RE->TC GABA_A decay (ms). See --sweep-tau2-re-tc.")
    ap.add_argument("--gh-tc", type=float, default=0.0,
                    help="TC's Ca2+-dependent I_h (ihca.mod) conductance. "
                         "Declared but never actually wired to TCCell until "
                         "now (silently ignored as gh=0.0). Candidate lever "
                         "for RE/TC multi-cycle waxing/waning: article "
                         "mechanism is Ca2+ from successive rebound bursts "
                         "locking I_h open, eventually terminating the "
                         "spindle and setting the ~2.5 s refractory period. "
                         "Default 0.0 keeps it off (untested territory).")
    a = ap.parse_args()

    if a.sweep_g_re_re:
        pc = h.ParallelContext()
        rank, nhost = int(pc.id()), int(pc.nhost())
        sizes = {k: max(1, int(round(v * a.scale))) for k, v in DEFAULT_SIZES.items()}
        if rank == 0:
            print("=== g_re_re sweep (%d ranks, scale=%.2f, tstop=%.0f ms) ==="
                  % (nhost, a.scale, a.sweep_tstop))
            print("%-10s %-9s %-10s %-12s %-10s %-10s"
                  % ("g_re_re", "RE spikes", "frac_burst", "mean_burst",
                     "n_events", "event_Hz"))
            print("-" * 62)
        for g_re_re in [float(x) for x in a.sweep_g_re_re.split(",")]:
            net = ParallelCorticoThalamicNet(sizes=sizes, conv=a.conv,
                                              het=a.het,
                                              delay_jitter=a.delay_jitter,
                                              g_re_re=g_re_re)
            net.run(tstop=a.sweep_tstop)
            t, g = net.gather()
            if rank == 0:
                re_lo, re_hi = net.ranges["re"]
                s = _re_burst_stats(t, g, re_lo, re_hi)
                print("%-10.4f %-9d %-10.3f %-12.2f %-10d %-10.3f"
                      % (g_re_re, s["n_spikes"], s["frac_burst"],
                         s["mean_burst_size"], s["n_events"], s["event_hz"]))
            net.teardown()
            del net
        if rank == 0:
            print("-" * 62)
        pc.barrier()
        pc.done()
        h.quit()
        return

    if a.sweep_g_l5_rec:
        pc = h.ParallelContext()
        rank, nhost = int(pc.id()), int(pc.nhost())
        sizes = {k: max(1, int(round(v * a.scale))) for k, v in DEFAULT_SIZES.items()}
        if rank == 0:
            print("=== g_l5_rec sweep (%d ranks, scale=%.2f, tstop=%.0f ms, "
                  "tau2_l5_rec=%.0f) ===" % (nhost, a.scale, a.sweep_tstop, a.tau2_l5_rec))
            print("%-9s %-9s %-8s %-9s %-9s | %-9s %-8s %-9s %-9s"
                  % ("g_l5_rec", "L5 spks", "L5frac_b", "L5burst", "L5hz",
                     "RE spks", "REfrac_b", "REburst", "REhz"))
            print("-" * 90)
        for g_l5_rec in [float(x) for x in a.sweep_g_l5_rec.split(",")]:
            net = ParallelCorticoThalamicNet(sizes=sizes, conv=a.conv,
                                              het=a.het,
                                              delay_jitter=a.delay_jitter,
                                              g_i_e_l5=a.g_i_e_l5,
                                              g_l5_rec=g_l5_rec,
                                              tau2_l5_rec=a.tau2_l5_rec,
                                              l5_rec_mech=a.l5_rec_mech,
                                              mg_l5_rec=a.mg_l5_rec)
            wall = net.run(tstop=a.sweep_tstop)
            t, g = net.gather()
            if rank == 0:
                l5_lo, l5_hi = net.ranges["l5e"]
                re_lo, re_hi = net.ranges["re"]
                sl5 = _re_burst_stats(t, g, l5_lo, l5_hi)
                sre = _re_burst_stats(t, g, re_lo, re_hi)
                print("%-9.4f %-9d %-8.3f %-9.2f %-9.3f | %-9d %-8.3f %-9.2f %-9.3f  (wall %.0fs)"
                      % (g_l5_rec, sl5["n_spikes"], sl5["frac_burst"],
                         sl5["mean_burst_size"], sl5["event_hz"],
                         sre["n_spikes"], sre["frac_burst"],
                         sre["mean_burst_size"], sre["event_hz"], wall))
            net.teardown()
            del net
        if rank == 0:
            print("-" * 90)
        pc.barrier()
        pc.done()
        h.quit()
        return

    if a.sweep_tau2_re_tc:
        pc = h.ParallelContext()
        rank, nhost = int(pc.id()), int(pc.nhost())
        sizes = {k: max(1, int(round(v * a.scale))) for k, v in DEFAULT_SIZES.items()}
        if rank == 0:
            print("=== tau2_re_tc sweep (%d ranks, scale=%.2f, tstop=%.0f ms) ==="
                  % (nhost, a.scale, a.sweep_tstop))
            print("%-9s %-9s %-8s %-9s %-9s | %-9s %-8s %-9s %-9s"
                  % ("tau2", "TC spks", "TCfrac_b", "TCburst", "TChz",
                     "RE spks", "REfrac_b", "REburst", "REhz"))
            print("-" * 90)
        for tau2 in [float(x) for x in a.sweep_tau2_re_tc.split(",")]:
            net = ParallelCorticoThalamicNet(sizes=sizes, conv=a.conv,
                                              het=a.het,
                                              delay_jitter=a.delay_jitter,
                                              g_i_e_l5=a.g_i_e_l5,
                                              tau2_re_tc=tau2)
            wall = net.run(tstop=a.sweep_tstop)
            t, g = net.gather()
            if rank == 0:
                tc_lo, tc_hi = net.ranges["tc"]
                re_lo, re_hi = net.ranges["re"]
                stc = _re_burst_stats(t, g, tc_lo, tc_hi)
                sre = _re_burst_stats(t, g, re_lo, re_hi)
                print("%-9.1f %-9d %-8.3f %-9.2f %-9.3f | %-9d %-8.3f %-9.2f %-9.3f  (wall %.0fs)"
                      % (tau2, stc["n_spikes"], stc["frac_burst"],
                         stc["mean_burst_size"], stc["event_hz"],
                         sre["n_spikes"], sre["frac_burst"],
                         sre["mean_burst_size"], sre["event_hz"], wall))
            net.teardown()
            del net
        if rank == 0:
            print("-" * 90)
        pc.barrier()
        pc.done()
        h.quit()
        return

    if a.sweep_gh_tc:
        pc = h.ParallelContext()
        rank, nhost = int(pc.id()), int(pc.nhost())
        sizes = {k: max(1, int(round(v * a.scale))) for k, v in DEFAULT_SIZES.items()}
        if rank == 0:
            print("=== gh_tc sweep (%d ranks, scale=%.2f, tstop=%.0f ms) ==="
                  % (nhost, a.scale, a.sweep_tstop))
            print("%-9s %-9s %-8s %-9s %-9s | %-9s %-8s %-9s %-9s"
                  % ("gh_tc", "TC spks", "TCfrac_b", "TCburst", "TChz",
                     "RE spks", "REfrac_b", "REburst", "REhz"))
            print("-" * 90)
        for gh_tc in [float(x) for x in a.sweep_gh_tc.split(",")]:
            net = ParallelCorticoThalamicNet(sizes=sizes, conv=a.conv,
                                              het=a.het,
                                              delay_jitter=a.delay_jitter,
                                              g_i_e_l5=a.g_i_e_l5, gh_tc=gh_tc)
            wall = net.run(tstop=a.sweep_tstop)
            t, g = net.gather()
            if rank == 0:
                tc_lo, tc_hi = net.ranges["tc"]
                re_lo, re_hi = net.ranges["re"]
                stc = _re_burst_stats(t, g, tc_lo, tc_hi)
                sre = _re_burst_stats(t, g, re_lo, re_hi)
                print("%-9.1e %-9d %-8.3f %-9.2f %-9.3f | %-9d %-8.3f %-9.2f %-9.3f  (wall %.0fs)"
                      % (gh_tc, stc["n_spikes"], stc["frac_burst"],
                         stc["mean_burst_size"], stc["event_hz"],
                         sre["n_spikes"], sre["frac_burst"],
                         sre["mean_burst_size"], sre["event_hz"], wall))
            net.teardown()
            del net
        if rank == 0:
            print("-" * 90)
        pc.barrier()
        pc.done()
        h.quit()
        return

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
                                      delay_jitter=a.delay_jitter,
                                      g_re_re=a.g_re_re, g_re_re_sd=a.g_re_re_sd,
                                      g_i_e_l5=a.g_i_e_l5, gh_tc=a.gh_tc,
                                      tau2_re_tc=a.tau2_re_tc,
                                      g_l5_rec=a.g_l5_rec, tau2_l5_rec=a.tau2_l5_rec,
                                      l5_rec_mech=a.l5_rec_mech, mg_l5_rec=a.mg_l5_rec)
    wall = net.run(tstop=a.tstop)
    t, g = net.gather()
    if net.rank == 0:
        _report(net, t, g, wall, a.tstop)
        re_lo, re_hi = net.ranges["re"]
        s = _re_burst_stats(t, g, re_lo, re_hi)
        print("RE burst shape (full run): frac_burst=%.3f mean_burst_size=%.2f "
              "n_events=%d event_Hz=%.3f (g_re_re=%.4f +/- %.4f)"
              % (s["frac_burst"], s["mean_burst_size"], s["n_events"],
                 s["event_hz"], a.g_re_re, a.g_re_re_sd))
        if a.tstop > 40000:
            # temporal stability check: does burst shape hold up over the
            # full run, or drift toward runaway/silence in the second half?
            mid = a.tstop / 2.0
            m1, m2 = t < mid, t >= mid
            s1 = _re_burst_stats(t[m1], g[m1], re_lo, re_hi)
            s2 = _re_burst_stats(t[m2], g[m2], re_lo, re_hi)
            print("  first half:  frac_burst=%.3f mean_burst_size=%.2f event_Hz=%.3f"
                  % (s1["frac_burst"], s1["mean_burst_size"], s1["event_hz"]))
            print("  second half: frac_burst=%.3f mean_burst_size=%.2f event_Hz=%.3f"
                  % (s2["frac_burst"], s2["mean_burst_size"], s2["event_hz"]))
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
