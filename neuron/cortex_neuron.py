"""
NEURON cortical column cell types: regular-spiking pyramidal cells (PYCell)
and fast-spiking PV+ interneurons (FSCell), plus CorticalColumn -- the L4 /
L2/3 / L5 / L6 laminar assembly wired per Mushtaq et al. 2024 (Table 3) /
config/network_auditory_hh.yaml, replacing the sinusoidal cortex proxy in
tc_network_nrn.py with a real conductance-based column.

Biophysics (Fernandez & Luthi 2020, sect. VI):
  PYCell   -- regular-spiking pyramidal: hh2 (fast Na+/K+ spikes) + ical (HVA
              Ca2+, opens on every spike upstroke) + cad (private submembrane
              Ca2+ pool) + sk2 (SK2 Ca2+-activated K+). This is the SAME SK2 +
              Ca2+-pool pair already validated in the thalamic RE cell
              (mod/sk2.mod, mod/cad.mod), reused here as the cortical
              spike-frequency adaptation current.
  PYCellIB -- intrinsically-bursting L5 pyramidal ("TuftIB", tc_architecture.py):
              PYCell + inap (persistent Na+, mod/inap.mod) and a SLOWED
              cad/sk2 recovery (taur 80ms -> 500ms). I_NaP turns a single
              spike into a self-sustaining plateau burst; the slow SK2/Ca2+
              recovery is what turns that from within-burst adaptation into a
              genuine ~1 Hz UP/DOWN relaxation oscillator. This is the
              column's OWN slow-oscillation generator -- no externally
              injected volley is needed for L5 to produce periodic UP states,
              replacing the imposed 1 Hz sinusoid tc_network_nrn.py used
              previously.
  FSCell   -- hh2 only, high gkbar (fast-spiking, minimal adaptation): the
              PV+ basket cell that recruits fast feedforward inhibition onto
              pyramidal somata (review sect. VI: TC spindle bursts recruit
              PV+ interneurons via GluA2-lacking AMPARs).

Layout (review sect. VI + Mushtaq Table 3, config/network_auditory_hh.yaml):
    TC (MGB) -> L4 (core, focal)          thalamocortical drive (sensory)
    L4  -> L2/3 -> L5 -> L6               intracortical excitatory flow
    L5  -> L2/3                           apical-tuft feedback (entrains
                                           superficial layers to the L5 SO)
    each layer: E -> I (feedforward) and I -> E (feedback inhibition)
    L6 -> TC, L6 -> RE                    corticothalamic feedback (closes
                                           the loop; the "L6 kick" that
                                           triggers spindles, review sect. V.B.2)
L5's PYCellIB population paces the column at ~1.4 Hz intrinsically; L4/L2-3
otherwise stay quiet without sensory/thalamic drive (see
ctx_thalamus_network.py for the TC -> L4 connection that activates them).

    ../.venv-neuron/bin/python cortex_neuron.py   # single-cell + column demo
"""

import os
import numpy as np
from neuron import h

h.load_file("stdrun.hoc")
_here = os.path.dirname(os.path.abspath(__file__))
if not hasattr(h, "ical"):
    _found = False
    for arch in ("arm64", "x86_64", "aarch64", "."):
        for lib in ("libnrnmech.so", "libnrnmech.dylib", ".libs/libnrnmech.so"):
            dll = os.path.join(_here, arch, lib)
            if os.path.exists(dll):
                h.nrn_load_dll(dll)
                _found = True
                break
        if _found:
            break
    if not _found and not hasattr(h, "ical"):
        raise RuntimeError(
            "NMODL mechanisms not found under %s.\nCompile them first:\n"
            "    cd %s && nrnivmodl mod" % (_here, _here))


class PYCell:
    """Regular-spiking cortical pyramidal cell.

    Two compartments (soma + apical dend) so thalamocortical/intracortical
    synapses can target the dendrite while spikes are generated at the soma,
    mirroring TCCell2C's electrotonic-separation design. gsk/gh default on:
    the SK2 + Ca2+ pool is the cortical adaptation current this model adds.
    """

    def __init__(self, gna=0.1, gk=0.012, gca=1.2e-4, gsk=1.5e-4,
                 depth=5.0, kd=5e-4, taur=80.0, e_pas=-68.0, dend_L=300.0,
                 dend_diam=3.0, Ra=150.0, gnap=0.0):
        self.soma = h.Section(name="pysoma", cell=self)
        self.soma.L = self.soma.diam = 25
        self.soma.Ra, self.soma.cm = 100, 1
        self.soma.insert("hh2")
        self.soma.gnabar_hh2, self.soma.gkbar_hh2 = gna, gk
        self.soma.insert("pas")
        self.soma.g_pas, self.soma.e_pas = 5e-5, e_pas
        # I_NaP: intrinsic bursting (L5 "TuftIB"; Compte et al. 2003 /
        # Bazhenov-Timofeev cortical SO models). Off by default (RS cell);
        # PYCellIB below turns it on. A persistent inward current that turns
        # a single somatic spike into a self-sustaining depolarising plateau
        # -- the cortical column's OWN slow-oscillation UP-state generator,
        # terminated by the SK2/cad adaptation already on the dendrite.
        if gnap > 0:
            self.soma.insert("inap")
            self.soma.gnabar_inap = gnap

        self.dend = h.Section(name="pydend", cell=self)
        self.dend.L, self.dend.diam = dend_L, dend_diam
        self.dend.Ra, self.dend.cm = Ra, 1
        self.dend.nseg = 5
        self.dend.insert("pas")
        self.dend.g_pas, self.dend.e_pas = 3e-5, e_pas
        if gca > 0:
            self.dend.insert("ical")
            self.dend.gcabar_ical = gca
            if gsk > 0:
                self.dend.insert("cad")
                self.dend.depth_cad, self.dend.taur_cad = depth, taur
                self.dend.insert("sk2")
                self.dend.gkbar_sk2, self.dend.kd_sk2 = gsk, kd
        self.dend.connect(self.soma(1))
        self.synsec = self.dend
        self._gaps = []

    def record(self):
        self.t = h.Vector().record(h._ref_t)
        self.vsoma = h.Vector().record(self.soma(0.5)._ref_v)
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        nc.threshold = -10
        self.spikes = h.Vector()
        nc.record(self.spikes)
        self._nc = nc


class PYCellIB(PYCell):
    """Intrinsically-bursting L5 pyramidal cell ("TuftIB", tc_architecture.py):
    PYCell + I_NaP (mod/inap.mod). A single suprathreshold event turns into a
    self-sustaining depolarising plateau (a burst of ~5 spikes), terminated by
    the SK2/Ca2+ adaptation on the dendrite. `taur` (the submembrane Ca2+
    pool's clearance time constant, mod/cad.mod) is slowed from the RS cell's
    80 ms to 500 ms, which is what turns SK2 from a fast within-burst
    adaptation current into a SLOW recovery variable -- the difference
    between spike-frequency adaptation and a genuine ~1 Hz UP/DOWN
    relaxation oscillator. This IS the cortical column's own
    slow-oscillation generator: no externally injected volley is needed for
    the column to produce periodic UP states."""

    def __init__(self, gnap=2e-4, gsk=8e-4, taur=500.0, **kwargs):
        super().__init__(gnap=gnap, gsk=gsk, taur=taur, **kwargs)


class FSCell:
    """Fast-spiking PV+ basket interneuron: hh2 only, high gkbar for narrow
    spikes and minimal adaptation (no Ca2+/SK2 -- PV+ cells are essentially
    non-adapting)."""

    def __init__(self, gna=0.1, gk=0.02, e_pas=-70.0):
        self.soma = h.Section(name="fs", cell=self)
        self.soma.L = self.soma.diam = 15
        self.soma.Ra, self.soma.cm = 100, 1
        self.soma.insert("hh2")
        self.soma.gnabar_hh2, self.soma.gkbar_hh2 = gna, gk
        self.soma.insert("pas")
        self.soma.g_pas, self.soma.e_pas = 5e-5, e_pas
        self.synsec = self.soma
        self._gaps = []

    def record(self):
        self.t = h.Vector().record(h._ref_t)
        self.vsoma = h.Vector().record(self.soma(0.5)._ref_v)
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        nc.threshold = -10
        self.spikes = h.Vector()
        nc.record(self.spikes)
        self._nc = nc


# Population sizes, config/network_auditory_hh.yaml ("Mushtaq et al. 2024",
# PY 200 / IN 40 total split across L4/L2/3/L5/L6). Scaled down 5x here for a
# fast NEURON demo; ratios preserved.
LAYER_SIZES = {
    "L4":  {"E": 10, "I": 2},
    "L23": {"E": 12, "I": 3},
    "L5":  {"E": 10, "I": 2},
    "L6":  {"E": 8,  "I": 1},
}


class CorticalColumn:
    """One auditory-cortex column: L4 / L2/3 / L5 / L6, each E (PYCell) + I
    (FSCell), wired per Mushtaq Table 3 / config/network_auditory_hh.yaml:
    L4 -> L2/3 -> L5 -> L6 (feedforward), L6 <-> L5 (recurrent), E -> I
    (feedforward) and I -> E (feedback inhibition) within each layer.
    """

    def __init__(self, sizes=None, seed=2, g_ff=0.0015, g_rec=0.0,
                 g_e_i=0.02, g_i_e=0.08, gsk=8e-4, e_pas=-70.0,
                 ib_frac=None):
        self.rng = np.random.default_rng(seed)
        sizes = sizes or LAYER_SIZES
        # Fraction of each layer's E population that is intrinsically
        # bursting (PYCellIB) rather than regular-spiking (PYCell).
        # tc_architecture.py labels L5 "TuftRS + TuftIB": L5 is where the
        # column's own slow oscillation originates (Steriade/Bazhenov-style
        # cortical SO models), so it gets the IB cells; L4/L2-3/L6 stay
        # purely regular-spiking.
        ib_frac = ib_frac or {"L5": 0.5}
        self.layers = {}
        for name, n in sizes.items():
            frac = ib_frac.get(name, 0.0)
            n_ib = int(round(n["E"] * frac))
            E = ([PYCellIB(e_pas=e_pas) for _ in range(n_ib)] +
                 [PYCell(e_pas=e_pas, gsk=gsk) for _ in range(n["E"] - n_ib)])
            I = [FSCell(e_pas=e_pas + 2) for _ in range(n["I"])]
            self.layers[name] = {"E": E, "I": I}
        self._syn, self._nc = [], []

        self._wire_intracortical(g_ff, g_rec)
        for name in self.layers:
            self._wire_layer_inh(name, g_e_i, g_i_e)

    # -- helpers ------------------------------------------------------------
    def _connect(self, pre, post, e, tau1, tau2, w, delay=1.0, k=None):
        tgt = getattr(post, "synsec", post.soma)
        syn = h.Exp2Syn(tgt(0.5)); syn.e = e; syn.tau1 = tau1; syn.tau2 = tau2
        nc = h.NetCon(pre.soma(0.5)._ref_v, syn, sec=pre.soma)
        nc.threshold = -10; nc.weight[0] = w if k is None else w / k
        nc.delay = delay
        self._syn.append(syn); self._nc.append(nc)

    def _project(self, pre_pop, post_pop, g, e=0, tau1=0.5, tau2=2.0, frac=0.5):
        """Sparse, fan-in-normalised excitatory projection pre_pop -> post_pop."""
        k = max(1, int(len(pre_pop) * frac))
        for post in post_pop:
            for j in self.rng.choice(len(pre_pop), k, replace=False):
                self._connect(pre_pop[j], post, e, tau1, tau2, g, k=k)

    # -- wiring ---------------------------------------------------------
    def _wire_intracortical(self, g_ff, g_rec):
        L4E, L23E = self.layers["L4"]["E"], self.layers["L23"]["E"]
        L5E, L6E = self.layers["L5"]["E"], self.layers["L6"]["E"]
        self._project(L4E, L23E, g_ff)     # L4 -> L2/3
        self._project(L23E, L5E, g_ff)     # L2/3 -> L5
        self._project(L4E, L5E, g_ff * 0.5)  # direct L4 -> L5
        # L5 -> L6 needs to be MUCH stronger than the other feedforward
        # links: L5's IB cells (PYCellIB) fire a brief, sparse burst (a
        # handful of cells, a few ms wide) once per SO cycle, not the
        # sustained volleys the other feedforward links carry, so the same
        # g_ff that works elsewhere leaves L6 below threshold most cycles.
        # At g_ff*6 every L5 burst reliably reaches L6 -- this is the "L6
        # kick" link the corticothalamic loop depends on.
        self._project(L5E, L6E, g_ff * 6.0)  # L5 -> L6
        if g_rec > 0:
            self._project(L6E, L5E, g_rec)   # L6 <-> L5 recurrent (optional;
            # feeding back into L5 at more than a token strength desynchronises
            # the IB cells' own SO clock into irregular double-bursting)
        # L5 -> L2/3 feedback (apical tuft projection): entrains the
        # superficial layers to the L5-generated slow oscillation even
        # without external/thalamic drive to L4.
        self._project(L5E, L23E, g_ff * 1.0)

    def _wire_layer_inh(self, name, g_e_i, g_i_e):
        E, I = self.layers[name]["E"], self.layers[name]["I"]
        self._project(E, I, g_e_i, e=0, tau1=0.5, tau2=2.0, frac=0.7)
        self._project(I, E, g_i_e, e=-75, tau1=1.0, tau2=8.0, frac=0.7)

    # -- introspection ----------------------------------------------------
    def all_cells(self):
        out = []
        for L in self.layers.values():
            out += L["E"] + L["I"]
        return out

    def record(self):
        for c in self.all_cells():
            c.record()


def _adaptation_demo():
    """Does a step current into a PYCell give spike-frequency adaptation via
    Ca2+ influx (ical) -> submembrane pool (cad) -> SK2?"""
    print("Cortical pyramidal spike-frequency adaptation (ical + cad + sk2):\n")
    print(f"{'gkbar_sk2':<12}{'n spikes':<10}{'ISI 1st (ms)':<14}{'ISI last (ms)':<14}")
    print("-" * 50)
    for gsk in [0.0, 5e-5, 1.5e-4, 4e-4]:
        cell = PYCell(gsk=gsk)
        cell.record()
        ic = h.IClamp(cell.soma(0.5)); ic.delay, ic.dur, ic.amp = 200, 800, 0.12
        h.celsius = 36
        h.finitialize(-68)
        h.continuerun(1200)
        sp = np.asarray(cell.spikes)
        sp = sp[(sp > ic.delay) & (sp < ic.delay + ic.dur)]
        isi = np.diff(sp)
        i0 = isi[0] if len(isi) else 0.0
        i1 = isi[-1] if len(isi) else 0.0
        print(f"{gsk:<12.1e}{len(sp):<10}{i0:<14.1f}{i1:<14.1f}")
    print("-" * 50)
    print("Rising gkbar_sk2 widens the late ISI relative to the first -- Ca2+-")
    print("dependent spike-frequency adaptation, the same SK2 + Ca2+ pool")
    print("mechanism validated in the thalamic RE cell, now in cortex.\n")


def _fs_demo():
    print("Fast-spiking PV+ interneuron (hh2 only, minimal adaptation):\n")
    cell = FSCell()
    cell.record()
    ic = h.IClamp(cell.soma(0.5)); ic.delay, ic.dur, ic.amp = 200, 800, 0.3
    h.celsius = 36
    h.finitialize(-70)
    h.continuerun(1200)
    sp = np.asarray(cell.spikes)
    sp = sp[(sp > ic.delay) & (sp < ic.delay + ic.dur)]
    isi = np.diff(sp)
    print(f"  {len(sp)} spikes, ISI ratio last/first = "
          f"{(isi[-1]/isi[0] if len(isi) else 0):.2f} (PV+ cells barely adapt)\n")


def _column_demo(tstop=8000.0):
    """L5's PYCellIB population (I_NaP + slow SK2/Ca2+ recovery) paces the
    column on its own -- no externally injected volley. Look for a periodic,
    self-generated UP state propagating L5 -> L6 (and L5 -> L2/3)."""
    print("Cortical column (L4/L2/3/L5/L6): self-generated slow oscillation,")
    print("no external drive -- L5's intrinsically-bursting cells are the")
    print("column's own UP-state pacemaker.\n")
    col = CorticalColumn()
    col.record()
    h.celsius = 36
    h.finitialize(-70)
    h.continuerun(tstop)
    for name in ["L4", "L23", "L5", "L6"]:
        sp = np.sort(np.concatenate([np.asarray(c.spikes) for c in col.layers[name]["E"]]))
        n = len(sp)
        rate = n / (len(col.layers[name]["E"]) * tstop / 1000.0)
        bursts = sp[np.diff(sp, prepend=-1e9) > 50] if n else np.array([])
        period = np.mean(np.diff(bursts)) if len(bursts) > 1 else 0.0
        print(f"  {name} E: {n} spikes, {rate:.1f} Hz/cell mean, "
              f"{len(bursts)} UP states, inter-UP-state interval ~{period:.0f} ms")
    print("\n(L4/L2/3 need thalamic/sensory drive to fire -- see")
    print("ctx_thalamus_network.py for TC -> L4.)")


if __name__ == "__main__":
    _adaptation_demo()
    _fs_demo()
    _column_demo()
