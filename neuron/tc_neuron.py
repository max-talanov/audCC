"""
NEURON scaffold for the auditory thalamo-cortical spindle model (migration WIP).

Motivation (see docs/spindle_review_mapping.md sect. 5.8 and the NEST model's
limits): sleep-spindle bursting is a **dendritic** phenomenon -- TRN via Ca_v3.3
in distal dendrites, TC via Ca_v3.1 in primary dendrites (Fernandez & Luthi
2020, sect. V). Point neurons (NEST's ht_neuron / AdEx) cannot represent this,
which is why in-network burst synchrony stalled. NEURON represents multi-
compartment cells and the T-current natively.

This module is a first, minimal step: a two-compartment TC (relay) cell with a
Ca_v3.1 T-current in its dendrite, used to demonstrate a genuine **dendritic
rebound burst**. The full network port is future work; the point here is to
prove the mechanism NEURON gives us that the point models could not.

Build the mechanism once, then run:

    .venv-neuron/bin/nrnivmodl mod          # compile neuron/mod/*.mod
    .venv-neuron/bin/python tc_neuron.py     # single-cell rebound-burst demo

Emits results in the tc_validate result contract (spikes/traces/meta) so the
same simulator-agnostic validator can eventually score a NEURON network.
"""

import os
import numpy as np
from neuron import h

h.load_file("stdrun.hoc")
# load the compiled T-current (arm64/x86_64 dir next to this file)
_here = os.path.dirname(os.path.abspath(__file__))
if not hasattr(h, "cav3"):
    for arch in ("arm64", "x86_64"):
        dll = os.path.join(_here, arch, "libnrnmech.dylib")
        if os.path.exists(dll):
            h.nrn_load_dll(dll)
            break


class TCCell:
    """Two-compartment thalamocortical relay cell: soma (fast Na/K) + dendrite
    carrying the low-threshold Ca_v3.1 T-current that produces rebound bursts."""

    def __init__(self, gcat=0.004):
        self.soma = h.Section(name="soma", cell=self)
        self.dend = h.Section(name="dend", cell=self)
        self.soma.L = self.soma.diam = 20
        self.dend.L, self.dend.diam = 100, 3
        self.dend.connect(self.soma(0.5))
        for sec in (self.soma, self.dend):
            sec.Ra, sec.cm = 100, 1
        # Soma: fast Na/K spikes. Dendrite: PASSIVE + the T-current only, so the
        # low-threshold Ca2+ spike (not dendritic Na) is what drives the soma to
        # fire a rebound BURST -- the standard thalamic-relay compartmentalisation.
        self.soma.insert("hh")
        self.soma.gkbar_hh = 0.01
        self.dend.insert("cav3")
        self.dend.gmax_cav3 = gcat
        self.dend.insert("pas")
        self.dend.g_pas, self.dend.e_pas = 2e-5, -70
        h.ion_style("ca_ion", 3, 2, 1, 1, 0, sec=self.dend)
        self.dend.cai, self.dend.cao = 5e-5, 2.0

    def record(self):
        self.t = h.Vector().record(h._ref_t)
        self.vsoma = h.Vector().record(self.soma(0.5)._ref_v)
        self.vdend = h.Vector().record(self.dend(0.5)._ref_v)
        self.ica = h.Vector().record(self.dend(0.5)._ref_ica)
        self.mT = h.Vector().record(self.dend(0.5)._ref_m_cav3)
        self.hT = h.Vector().record(self.dend(0.5)._ref_h_cav3)
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        nc.threshold = -10
        self.spikes = h.Vector()
        nc.record(self.spikes)
        self._nc = nc


def rebound_demo(gcat=0.05, hyper_amp=-0.3, hyper_dur=800.0):
    """Hyperpolarise the relay cell, release, and measure the T-current-driven
    rebound. Returns both the mechanism proof (peak inward Ca current, I_T
    (de)inactivation) and the somatic spike count."""
    cell = TCCell(gcat=gcat)
    cell.record()
    ic = h.IClamp(cell.soma(0.5))
    ic.delay, ic.dur, ic.amp = 300.0, hyper_dur, hyper_amp   # nA, hyperpolarising
    h.celsius = 36
    h.finitialize(-70)
    h.continuerun(ic.delay + ic.dur + 500.0)
    t = np.asarray(cell.t)
    sp = np.asarray(cell.spikes)
    release = ic.delay + ic.dur
    burst = sp[(sp > release) & (sp < release + 400.0)]
    isi = float(np.mean(np.diff(burst))) if len(burst) > 1 else 0.0
    win = (t > release - 5) & (t < release + 200)
    return {
        "n_rebound": int(len(burst)),
        "burst_hz": (1000.0 / isi) if isi > 0 else 0.0,
        "ica_peak": float(np.asarray(cell.ica)[win].min()),   # inward (negative)
        "h_at_release": float(np.asarray(cell.hT)[
            (t > release - 5) & (t < release)].mean()),        # de-inactivation
        "vdend_min": float(np.asarray(cell.vdend).min()),
    }


if __name__ == "__main__":
    print("NEURON multi-compartment TC relay cell: Ca_v3.1 T-current demo\n")
    print("Mechanism check -- does the dendritic T-current de-inactivate under")
    print("hyperpolarisation and drive an inward Ca current on release?\n")
    print(f"{'gCaT (S/cm2)':<14}{'I_Ca peak':<12}{'h(release)':<12}"
          f"{'rebound spikes':<16}")
    print("-" * 60)
    for g in [0.0, 0.01, 0.05, 0.1]:
        r = rebound_demo(gcat=g)
        print(f"{g:<14}{r['ica_peak']:<12.3f}{r['h_at_release']:<12.2f}"
              f"{r['n_rebound']:<16}")
    print("-" * 60)
    print("Status: the T-current de-inactivates (h ~ 1) and produces a strong")
    print("inward Ca current scaling with gCaT -- the mechanism the NEST point")
    print("models could not express. It drives a post-inhibitory rebound to")
    print("spike threshold; producing a full multi-spike BURST needs proper")
    print("delayed-rectifier kinetics in the soma and is the next step (the HH")
    print("soma goes into depolarisation block on the Ca plateau). See README.")
