"""
Full corticothalamic network: the cortical column (cortex_neuron.py) wired to
the thalamic TC<->RE loop (tc_neuron.py / tc_network_nrn.py), closing the
sleep-spindle loop with a REAL spiking cortex instead of the sinusoidal SO
proxy `tc_network_nrn.ThalamicNet._wire_cortex` used.

Wiring (Fernandez & Luthi 2020, sect. V-VI; Mushtaq et al. 2024 Table 3):
  L5 IB cells (PYCellIB) -> L6         cortex's OWN ~1.4 Hz slow-oscillation
                                        pacemaker (cortex_neuron.py) -- no
                                        externally injected volley
  L6 E -> TC, L6 E -> RE              corticothalamic feedback: this is the
                                       "L6 kick" that triggers spindles
                                       (review sect. V.B.2) -- genuine L6
                                       pyramidal spikes, not an injected
                                       NetStim clock.
  RE -> TC   GABA_A (+ optional GABA_B)
  TC -> RE   AMPA
  RE -> RE   GABA_A lateral inhibition + gap junctions
  TC -> L4 E  (optional, off by default) sensory/auditory drive

SK2 + Ca2+-dependent channels appear in BOTH structures, as in the article:
  thalamic RE  -- SK2 burst terminator (mod/sk2.mod + mod/cad.mod)
  thalamic TC  -- Ca2+-dependent I_h spindle terminator (mod/ihca.mod)
  cortical PY  -- SK2 spike-frequency adaptation (mod/sk2.mod + mod/cad.mod,
                  fed by the HVA Ca2+ current mod/ical.mod) -- the same
                  building block reused for the cortical UP-state / adaptation
                  mechanism (review sect. VI).

    ../.venv-neuron/bin/python ctx_thalamus_network.py   # full network demo
"""

import numpy as np
from neuron import h

h.load_file("stdrun.hoc")

try:
    from . import tc_neuron as T
    from . import cortex_neuron as C
except ImportError:
    import tc_neuron as T
    import cortex_neuron as C


class CorticoThalamicNet:
    def __init__(self, n_tc=10, n_re=10, seed=1,
                 g_re_tc=0.015, g_re_tc_b=0.0, g_tc_re=0.011,
                 g_re_re=0.002, g_gap=0.03,
                 tc_bias=0.0, re_bias=0.0, tc_e_pas=-80.0, re_e_pas=-82.0,
                 gsk_re=5e-5, gh_tc=0.0,
                 so_freq=1.0, g_l4_drive=0.03,
                 g_tc_l4=0.0005, g_l6_tc=0.03, g_l6_re=0.03,
                 het=0.05, delay_jitter=0.0,
                 column_kwargs=None):
        self.rng = np.random.default_rng(seed)
        self.n_tc, self.n_re = n_tc, n_re
        self.so_freq, self.g_l4_drive = so_freq, g_l4_drive
        # het: per-cell relative jitter (uniform +/- het) on the intrinsic
        # conductance densities that set burst threshold/timing (gcabar_it(2),
        # e_pas). delay_jitter: +/- ms uniform jitter on every synaptic delay.
        # Both are OFF by default (het=0) so nothing here changes prior
        # behaviour; see neuron/README.md "Heterogeneity" for what this tests.
        self.het, self.delay_jitter = het, delay_jitter

        # -- thalamus: same cells/calibration as tc_network_nrn.ThalamicNet --
        self.tc = [T.TCCell(gsk=0.0, gh=gh_tc) for _ in range(n_tc)]
        for c in self.tc:
            c.soma.e_pas = tc_e_pas * self._jitter(1.0)
            c.soma.gcabar_it *= self._jitter(1.0)
        self.re = [T.RECell(gsk=gsk_re) for _ in range(n_re)]
        for c in self.re:
            c.soma.e_pas = re_e_pas * self._jitter(1.0)
            c.soma.gcabar_it2 *= self._jitter(1.0)

        # -- cortex: the real column, replacing the sinusoidal SO proxy --
        self.column = C.CorticalColumn(seed=seed + 1, het=het,
                                        **(column_kwargs or {}))

        self._syn, self._nc, self._stim, self._gaps, self._ic = [], [], [], [], []

        for c in self.tc:
            self._bias(c, tc_bias)
        for c in self.re:
            self._bias(c, re_bias)

        self._wire_re_tc(g_re_tc, g_re_tc_b)
        self._wire_tc_re(g_tc_re)
        self._wire_re_re(g_re_re)
        self._wire_gap(g_gap)
        self._wire_thalamocortical(g_tc_l4)   # TC -> L4
        self._wire_corticothalamic(g_l6_tc, g_l6_re)  # L6 -> TC, L6 -> RE

    # -- helpers --------------------------------------------------------
    def _jitter(self, base):
        if self.het == 0:
            return base
        return base * (1.0 + self.rng.uniform(-self.het, self.het))

    def _bias(self, cell, amp):
        if amp == 0:
            return
        ic = h.IClamp(cell.soma(0.5)); ic.delay = 0; ic.dur = 1e9; ic.amp = amp
        self._ic.append(ic)

    def _syn_connect(self, pre, post, e, tau1, tau2, w, delay):
        tgt = getattr(post, "synsec", post.soma)
        syn = h.Exp2Syn(tgt(0.5)); syn.e = e; syn.tau1 = tau1; syn.tau2 = tau2
        nc = h.NetCon(pre.soma(0.5)._ref_v, syn, sec=pre.soma)
        nc.threshold = -10; nc.weight[0] = w
        if self.delay_jitter > 0:
            delay = max(0.1, delay + self.rng.uniform(-self.delay_jitter,
                                                        self.delay_jitter))
        nc.delay = delay
        self._syn.append(syn); self._nc.append(nc)

    def _project(self, pre_pop, post_pop, g, e, tau1, tau2, frac=0.5, delay=1.0):
        k = max(1, int(len(pre_pop) * frac))
        for post in post_pop:
            for j in self.rng.choice(len(pre_pop), k, replace=False):
                self._syn_connect(pre_pop[j], post, e, tau1, tau2, g / k, delay)

    # -- intrathalamic wiring (unchanged from tc_network_nrn.ThalamicNet) --
    def _wire_re_tc(self, g, g_b=0.0):
        k = max(1, self.n_re // 2)
        for tc in self.tc:
            for j in self.rng.choice(self.n_re, k, replace=False):
                self._syn_connect(self.re[j], tc, e=-85, tau1=1, tau2=8,
                                  w=g / k, delay=1.0)
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

    # -- thalamocortical loop --------------------------------------------
    def _wire_thalamocortical(self, g):
        """TC -> L4 E: the core, first-order thalamocortical projection."""
        L4E = self.column.layers["L4"]["E"]
        self._project(self.tc, L4E, g, e=0, tau1=0.5, tau2=2, frac=0.6)

    def _wire_corticothalamic(self, g_tc, g_re):
        """L6 E -> TC and L6 E -> RE: the 'L6 kick' that triggers spindles
        (review sect. V.B.2), now genuine L6 pyramidal spikes. L6 fires in
        brief ~10 ms UP-state bursts (cortex_neuron.CorticalColumn), so a
        modest, sparse fan-in (frac=0.3) delivers a phasic kick per SO cycle
        rather than a tonic barrage."""
        L6E = self.column.layers["L6"]["E"]
        self._project(L6E, self.tc, g_tc, e=0, tau1=0.5, tau2=2, frac=0.8)
        self._project(L6E, self.re, g_re, e=0, tau1=0.5, tau2=2, frac=0.8)

    # -- run --------------------------------------------------------------
    def record(self):
        for c in self.tc + self.re + self.column.all_cells():
            c.record()

    def run(self, tstop=8000.0, drive_L4=False):
        self.record()
        if drive_L4:
            # OPTIONAL sensory-like volley to L4 (auditory input), on top of
            # the column's own SO. Not needed for the SO/spindle loop itself:
            # L5's PYCellIB population (cortex_neuron.py) paces the column
            # intrinsically at ~1.4 Hz and drives L6 -> TC/RE every cycle
            # without any external trigger.
            period = 1000.0 / self.so_freq
            for t0 in np.arange(300.0, tstop, period):
                ns = h.NetStim(); ns.interval = 4; ns.number = 6; ns.start = t0
                ns.noise = 0.2
                self._stim.append(ns)
                for e in self.column.layers["L4"]["E"]:
                    syn = h.Exp2Syn(e.dend(0.5)); syn.e = 0; syn.tau1 = 0.5; syn.tau2 = 2
                    nc = h.NetCon(ns, syn); nc.weight[0] = self.g_l4_drive; nc.delay = 1
                    self._syn.append(syn); self._nc.append(nc)
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
            if not cells:
                return np.zeros_like(grid)
            vs = [np.interp(grid, np.asarray(c.t), np.asarray(c.vsoma))
                  for c in cells]
            return np.mean(vs, axis=0)

        spikes = {"MGB": pop(self.tc, 0), "nRT": pop(self.re, self.n_tc)}
        offset = self.n_tc + self.n_re
        for name, L in self.column.layers.items():
            spikes[name] = pop(L["E"] + L["I"], offset)
            offset += len(L["E"]) + len(L["I"])

        traces = {
            "MGB": {"time": grid, "voltage": meanvm(self.tc)},
            "nRT": {"time": grid, "voltage": meanvm(self.re)},
        }
        for name, L in self.column.layers.items():
            traces[name] = {"time": grid, "voltage": meanvm(L["E"])}
        meta = {"tstop": tstop, "seed": 1}
        return spikes, traces, meta


if __name__ == "__main__":
    print("NEURON corticothalamic network: cortical column (L4/L2/3/L5/L6,")
    print("L5 self-paced) + thalamic TC<->RE loop, closed via L6->TC/RE.\n")
    net = CorticoThalamicNet()
    spikes, traces, meta = net.run(tstop=8000.0)
    for lab, key in [("MGB (TC)", "MGB"), ("nRT (RE)", "nRT"),
                      ("L4", "L4"), ("L2/3", "L23"), ("L5", "L5"), ("L6", "L6")]:
        n = len(spikes[key]["times"])
        print(f"  {lab:<10}: {n} spikes ({n/(meta['tstop']/1000.0):.1f} Hz total)")

    def _events(times, senders, cell=None, gap=50.0):
        t = np.sort(times[senders == cell]) if cell is not None else np.sort(times)
        if len(t) < 2:
            return t, 0.0
        starts = t[np.diff(t, prepend=-1e9) > gap]
        period = np.mean(np.diff(starts)) if len(starts) > 1 else 0.0
        return starts, period

    re_starts, re_period = _events(spikes["nRT"]["times"], spikes["nRT"]["senders"],
                                    cell=net.n_tc)
    print(f"\n  RE burst events (cell 0): {len(re_starts)}, "
          f"inter-event interval ~{re_period:.0f} ms "
          f"({1000.0/re_period if re_period else 0:.2f} Hz event rate)")
    print("  Each RE burst is a single fast rebound (T-current) volley locked")
    print("  to the L6 corticothalamic kick -- NOT yet a sustained 10-15 Hz")
    print("  multi-cycle spindle train. See neuron/README.md.")
