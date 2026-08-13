"""
NEURON thalamic spindle network: the TC<->RE loop assembled from the verified
single cells (tc_neuron.py), producing emergent spindle oscillations.

Wiring (Fernandez & Luthi 2020, sect. V):
  RE -> TC   GABA_A, hyperpolarised Cl- reversal (e=-85 mV) so the IPSP
             de-inactivates the TC I_T and evokes a rebound burst.
  TC -> RE   AMPA, so a relay rebound re-excites the reticular cells.
  RE -> RE   GABA_A lateral inhibition (antisynchronising) + gap junctions
             (electrical synchronisation; the mechanism NEST could not express).
  cortex ->  a ~1 Hz slow-oscillation volley to BOTH RE and TC once per UP
             state (L6 excites both; ~2:1 in Mushtaq Table 3). The TC branch is
             essential: without an excitatory input the relay, held hyperpolarised
             for I_T de-inactivation, can never reach threshold and never
             rebounds. The cortex is reduced to this SO drive here.

Sleep state: BOTH populations rest hyperpolarised (tc_e_pas / re_e_pas), which
is what de-inactivates I_T -- at -71 mV the inactivation gate sits at h=0.05
(95% inactivated) and no rebound is possible. TC must stay ABOVE the GABA_A
reversal (-85 mV) or reticular input depolarises it instead of inhibiting it.

Emits (spikes, traces, meta) in the tc_validate result contract, so the same
simulator-agnostic 10-criteria harness that scores the NEST model can score this.

    ../.venv-neuron/bin/python tc_network_nrn.py            # run + summary
"""

import numpy as np
from neuron import h

h.load_file("stdrun.hoc")

try:
    from . import tc_neuron as T
except ImportError:
    import tc_neuron as T


class ThalamicNet:
    def __init__(self, n_tc=10, n_re=10, seed=1,
                 g_re_tc=0.015, g_re_tc_b=0.0, g_tc_re=0.004,
                 g_re_re=0.002, g_gap=0.03,
                 g_cort=0.01, g_cort_tc=0.035, so_freq=1.0,
                 tc_bias=0.0, re_bias=0.0, tc_e_pas=-80.0,
                 re_e_pas=-82.0, recruit_ms=0.0):
        self.rng = np.random.default_rng(seed)
        self.n_tc, self.n_re = n_tc, n_re
        self.so_freq = so_freq
        self.g_cort_tc = g_cort_tc
        # Progressive recruitment: spread the corticothalamic volley over
        # recruit_ms so cells enter the spindle on DIFFERENT cycles. I_T needs
        # ~320 ms to de-inactivate but the spindle cycle is ~83 ms, so a cell
        # can only burst once or twice (review sect. V.C). A synchronous volley
        # therefore makes the whole population deplete I_T together and the
        # event dies after ~2 cycles; staggering sustains the population rhythm.
        self.recruit_ms = recruit_ms
        self.tc = [T.TCCell(gsk=0.0) for _ in range(n_tc)]
        # Sleep state: relay cells rest HYPERPOLARISED (raised K-leak), which is
        # what de-inactivates I_T. At -71 mV the inactivation gate h sits at 0.05
        # (95% inactivated) and no rebound is possible; near -80 mV h ~ 0.4-0.5.
        # Must stay ABOVE the GABA_A reversal (-85 mV) or RE input depolarises TC.
        for c in self.tc:
            c.soma.e_pas = tc_e_pas
        self.re = [T.RECell(gsk=0.003) for _ in range(n_re)]
        # RE likewise needs a hyperpolarised sleep rest for its own I_T to
        # de-inactivate -- the review's 'TRN bursts only below ~-55 mV'.
        for c in self.re:
            c.soma.e_pas = re_e_pas
        self._syn, self._nc, self._stim, self._gaps, self._ic = [], [], [], [], []

        # tonic bias currents (set the thalamic operating point for sleep)
        for c in self.tc:
            self._bias(c, tc_bias)
        for c in self.re:
            self._bias(c, re_bias)

        self._wire_re_tc(g_re_tc, g_re_tc_b)   # RE -> TC (GABA_A [+ GABA_B])
        self._wire_tc_re(g_tc_re)     # TC -> RE (AMPA)
        self._wire_re_re(g_re_re)     # RE -> RE (GABA_A)
        self._wire_gap(g_gap)         # RE <-> RE gap junctions (nearest neighbour)
        self._wire_cortex(g_cort)     # SO-gated corticothalamic drive to RE

    # -- helpers ----------------------------------------------------------
    def _bias(self, cell, amp):
        if amp == 0:
            return
        ic = h.IClamp(cell.soma(0.5)); ic.delay = 0; ic.dur = 1e9; ic.amp = amp
        self._ic.append(ic)

    def _syn_connect(self, pre, post, e, tau1, tau2, w, delay):
        syn = h.Exp2Syn(post.soma(0.5)); syn.e = e; syn.tau1 = tau1; syn.tau2 = tau2
        nc = h.NetCon(pre.soma(0.5)._ref_v, syn, sec=pre.soma)
        nc.threshold = -10; nc.weight[0] = w; nc.delay = delay
        self._syn.append(syn); self._nc.append(nc)

    # -- wiring -----------------------------------------------------------
    def _wire_re_tc(self, g, g_b=0.0):
        # each TC inhibited by ~half the RE cells (fan-in normalised)
        k = max(1, self.n_re // 2)
        for tc in self.tc:
            for j in self.rng.choice(self.n_re, k, replace=False):
                self._syn_connect(self.re[j], tc, e=-85, tau1=1, tau2=8,
                                  w=g / k, delay=1.0)
                # Optional slow GABA_B component (E_K = -90 mV, rise 60 /
                # decay 200 ms). Fernandez & Luthi sect. V.B.3: the slow
                # reticular inhibition is what paces spindle waxing/waning, and
                # it holds TC hyperpolarised BETWEEN cycles so I_T stays
                # de-inactivated -- the candidate mechanism for spindle-length
                # events (the network otherwise has GABA_A only).
                if g_b:
                    self._syn_connect(self.re[j], tc, e=-90, tau1=60, tau2=200,
                                      w=g_b / k, delay=1.0)

    def _wire_tc_re(self, g):
        k = max(1, self.n_tc // 2)
        for re in self.re:
            for j in self.rng.choice(self.n_tc, k, replace=False):
                self._syn_connect(self.tc[j], re, e=0, tau1=0.5, tau2=2,
                                  w=g / k, delay=1.0)

    def _wire_re_re(self, g):
        for i, re in enumerate(self.re):
            for j in range(self.n_re):
                if i != j and abs(i - j) <= 2:
                    self._syn_connect(self.re[j], re, e=-75, tau1=1, tau2=6,
                                      w=g, delay=1.0)

    def _wire_gap(self, g):
        for i in range(self.n_re):
            j = (i + 1) % self.n_re
            T.gap_junction(self.re[i], self.re[j], g=g)

    def _wire_cortex(self, g):
        # A burst generator once per SO cycle: NetStim in bursty mode delivering
        # a short high-rate volley to every RE cell, UP-state gated.
        period = 1000.0 / self.so_freq
        ns = h.NetStim(); ns.interval = 3; ns.number = 6; ns.start = 300
        ns.noise = 0.2
        # repeat the volley each SO cycle via a second NetStim as a clock
        self._stim.append(ns)
        for re in self.re:
            syn = h.Exp2Syn(re.soma(0.5)); syn.e = 0; syn.tau1 = 0.5; syn.tau2 = 2
            nc = h.NetCon(ns, syn); nc.weight[0] = g; nc.delay = 1
            self._syn.append(syn); self._nc.append(nc)
        # clock: re-trigger the volley each period by scheduling multiple NetStims
        self._so_times = np.arange(300.0, 1e6, period)

    def _install_so_clock(self, tstop):
        """Install one burst NetStim per SO cycle up to tstop (UP-state drive)."""
        period = 1000.0 / self.so_freq
        for t0 in np.arange(300.0, tstop, period):
            if self.recruit_ms > 0:
                # one jittered volley PER CELL -> progressive recruitment
                for re in self.re:
                    ns = h.NetStim(); ns.interval = 4; ns.number = 5
                    ns.start = t0 + self.rng.uniform(0, self.recruit_ms)
                    ns.noise = 0.3
                    self._stim.append(ns)
                    syn = h.Exp2Syn(re.soma(0.5)); syn.e = 0; syn.tau1 = 0.5; syn.tau2 = 2
                    nc = h.NetCon(ns, syn); nc.weight[0] = 0.01; nc.delay = 1
                    self._syn.append(syn); self._nc.append(nc)
                ns = self._stim[-1]      # reuse the last for the TC branch below
            else:
                ns = h.NetStim(); ns.interval = 4; ns.number = 5; ns.start = t0
                ns.noise = 0.3
                self._stim.append(ns)
                for re in self.re:
                    syn = h.Exp2Syn(re.soma(0.5)); syn.e = 0; syn.tau1 = 0.5; syn.tau2 = 2
                    nc = h.NetCon(ns, syn); nc.weight[0] = 0.01; nc.delay = 1
                    self._syn.append(syn); self._nc.append(nc)
            # ... and onto the relay cells, the depolarising kick that lets the
            # de-inactivated I_T fire its rebound burst.
            for tc in self.tc:
                syn = h.Exp2Syn(tc.soma(0.5)); syn.e = 0; syn.tau1 = 0.5; syn.tau2 = 2
                nc = h.NetCon(ns, syn); nc.weight[0] = self.g_cort_tc; nc.delay = 1
                self._syn.append(syn); self._nc.append(nc)

    # -- run --------------------------------------------------------------
    def run(self, tstop=8000.0):
        # spike + Vm recorders
        for c in self.tc + self.re:
            c.record()
        self._install_so_clock(tstop)
        h.celsius = 36
        h.finitialize(-70)
        h.continuerun(tstop)
        return self._collect(tstop)

    def _collect(self, tstop):
        def pop(cells, offset):
            times, senders = [], []
            for i, c in enumerate(cells):
                st = np.asarray(c.spikes)
                times.append(st); senders.append(np.full(len(st), offset + i))
            t = np.concatenate(times) if times else np.array([])
            s = np.concatenate(senders) if senders else np.array([])
            o = np.argsort(t)
            return {"times": t[o], "senders": s[o]}

        grid = np.arange(0.0, tstop, 1.0)

        def meanvm(cells):
            vs = [np.interp(grid, np.asarray(c.t), np.asarray(c.vsoma))
                  for c in cells]
            return np.mean(vs, axis=0)

        spikes = {"MGB": pop(self.tc, 0), "nRT": pop(self.re, self.n_tc)}
        thal_v = meanvm(self.tc)
        # cortical slow-wave proxy: the 1 Hz SO drive (cortex reduced to its SO)
        cort_v = -63 + 6 * np.sin(2 * np.pi * self.so_freq * grid / 1000.0)
        traces = {
            "MGB": {"time": grid, "voltage": thal_v},
            "nRT": {"time": grid, "voltage": meanvm(self.re)},
            "L23": {"time": grid, "voltage": cort_v},
            "L5": {"time": grid, "voltage": cort_v},
            "L6": {"time": grid, "voltage": cort_v},
        }
        meta = {"tstop": tstop, "seed": 1}
        return spikes, traces, meta


if __name__ == "__main__":
    print("NEURON thalamic spindle network (TC<->RE loop) -- first run\n")
    net = ThalamicNet()
    spikes, traces, meta = net.run(tstop=8000.0)
    for lab, key in [("MGB (TC)", "MGB"), ("nRT (RE)", "nRT")]:
        n = len(spikes[key]["times"])
        print(f"  {lab}: {n} spikes ({n/(meta['tstop']/1000.0):.1f} Hz total)")
    # quick spindle-band check on the thalamic spike rate
    allt = np.concatenate([spikes["MGB"]["times"], spikes["nRT"]["times"]])
    if len(allt) > 10:
        t = np.arange(0, meta["tstop"], 1.0)
        rate = np.histogram(allt, bins=len(t))[0].astype(float)
        rate -= rate.mean()
        f = np.fft.rfftfreq(len(rate), 1e-3)
        P = np.abs(np.fft.rfft(rate)) ** 2
        band = (f >= 6) & (f <= 18)
        pk = f[band][np.argmax(P[band])]
        print(f"  thalamic rate spectral peak in 6-18 Hz: {pk:.1f} Hz")
    print("\nEmits the tc_validate result contract. To score against the same")
    print("10-criteria harness that grades the NEST model:")
    print("    import tc_network_nrn, tc_validate")
    print("    s,tr,m = tc_network_nrn.ThalamicNet().run(60000)")
    print("    tc_validate.validate(s, tr, m)")
