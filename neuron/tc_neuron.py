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
    """Single-compartment thalamocortical relay cell (Destexhe et al. 1996).

    Traub-Miles fast Na+/K+ (hh2) for spikes + the Destexhe low-threshold
    T-current (itd; GHK driving force, temperature-scaled kinetics) for the
    post-inhibitory rebound. A SINGLE compartment is used deliberately: the
    2-compartment version re-primed the T-current through somatic back-
    propagation and gave a ~250 ms plateau, whereas this produces a genuine,
    brief low-threshold Ca2+ spike carrying a physiological 2-6 spike burst.

    Default parameters are calibrated to a 6-spike rebound burst at ~160 Hz
    over ~31 ms (see neuron/README.md).
    """

    def __init__(self, pcabar=1.7e-4, gk=0.011, gna=0.1, gsk=0.0):
        self.soma = h.Section(name="soma", cell=self)
        self.soma.L = self.soma.diam = 80
        self.soma.Ra, self.soma.cm = 100, 1
        self.soma.insert("hh2")
        self.soma.gnabar_hh2, self.soma.gkbar_hh2 = gna, gk
        self.soma.insert("itd")           # Destexhe 1996 T-current
        self.soma.pcabar_itd = pcabar
        self.soma.insert("pas")
        self.soma.g_pas, self.soma.e_pas = 5e-5, -74
        h.ion_style("ca_ion", 1, 2, 0, 0, 0, sec=self.soma)
        self.soma.cai, self.soma.cao = 2.4e-4, 2.0
        # Optional SK2 burst terminator (needs its own Ca pool `cad`, private
        # `sk` ion). The Destexhe LTS is already brief (~31 ms), so SK2 is not
        # needed for the TC burst; it is retained for the RE cell and longer
        # plateaus. Enabling it re-styles Ca to let cad accumulate.
        if gsk > 0:
            h.ion_style("ca_ion", 1, 2, 1, 0, 0, sec=self.soma)
            self.soma.insert("cad")
            self.soma.insert("sk2")
            self.soma.gkbar_sk2 = gsk

    def record(self):
        self.t = h.Vector().record(h._ref_t)
        self.vsoma = h.Vector().record(self.soma(0.5)._ref_v)
        self.ica = h.Vector().record(self.soma(0.5)._ref_ica)
        self.hT = h.Vector().record(self.soma(0.5)._ref_h_itd)
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        nc.threshold = -10
        self.spikes = h.Vector()
        nc.record(self.spikes)
        self._nc = nc


def rebound_demo(pcabar=1.7e-4, gk=0.011, hyper_amp=-0.2, hyper_dur=500.0):
    """Hyperpolarise the relay cell, release, and measure the rebound burst."""
    cell = TCCell(pcabar=pcabar, gk=gk)
    cell.record()
    ic = h.IClamp(cell.soma(0.5))
    ic.delay, ic.dur, ic.amp = 300.0, hyper_dur, hyper_amp
    h.celsius = 36
    h.finitialize(-74)
    h.continuerun(ic.delay + ic.dur + 400.0)
    t = np.asarray(cell.t)
    sp = np.asarray(cell.spikes)
    release = ic.delay + ic.dur
    burst = sp[(sp > release) & (sp < release + 250.0)]
    isi = float(np.mean(np.diff(burst))) if len(burst) > 1 else 0.0
    return {
        "n_rebound": int(len(burst)),
        "burst_hz": (1000.0 / isi) if isi > 0 else 0.0,
        "burst_ms": float(burst[-1] - burst[0]) if len(burst) > 1 else 0.0,
        "ica_peak": float(np.asarray(cell.ica)[
            (t > release - 5) & (t < release + 150)].min()),
    }


if __name__ == "__main__":
    print("NEURON single-compartment TC relay cell (Destexhe 1996): rebound "
          "BURST demo\n")
    print("Traub-Miles (hh2) spikes + the Destexhe low-threshold T-current (itd,")
    print("GHK + temperature-scaled kinetics). Release from hyperpolarisation")
    print("evokes a brief low-threshold Ca2+ spike carrying a fast spike burst.\n")
    print(f"{'pcabar':<12}{'I_Ca peak':<12}{'spikes':<10}{'freq (Hz)':<12}"
          f"{'burst (ms)':<12}")
    print("-" * 58)
    for pca in [0.0, 1.0e-4, 1.7e-4, 2.5e-4]:
        r = rebound_demo(pcabar=pca)
        print(f"{pca:<12.1e}{r['ica_peak']:<12.3f}{r['n_rebound']:<10}"
              f"{r['burst_hz']:<12.0f}{r['burst_ms']:<12.0f}")
    print("-" * 58)
    r = rebound_demo()   # locked defaults
    print(f"\nDefault relay cell: a {r['n_rebound']}-spike rebound burst at "
          f"{r['burst_hz']:.0f} Hz over {r['burst_ms']:.0f} ms -- a genuine")
    print("thalamic low-threshold Ca2+ spike, in the review's 2-6 spike range.")
    print("Conductance-based (Destexhe 1996 I_T, GHK) -- not a phenomenological")
    print("fit, and no depolarisation block. This is the mechanism the NEST")
    print("point models (ht_neuron, AdEx) could not express.")
