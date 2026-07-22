"""
Closed-loop auditory thalamo-cortical column with sleep rhythms (NEST).

Architecture is adapted from the Cortical Column model
(https://github.com/max-talanov/cc/tree/main/optimization_nest,
``core/simulator.py``): a five-stage column of ``iaf_cond_exp`` neurons --
thalamus (TCR excitatory / nRT inhibitory), L4, L2/3 (RS + FRB), L5
(TuftRS + TuftIB), L6, with a shared L5/6 interneuron pool (Basket / LTS /
Axoaxonic) -- wired with ``pairwise_bernoulli`` connections and recorded with
spike recorders and multimeters.

This module is *not* an optimisation harness (cc's ``run_opt.py`` /
``methods/`` are dropped); it builds one column and runs it. Two biological
additions turn cc's feed-forward column into a sleeping auditory loop:

1. **A closed thalamo-cortical loop.** cc wires thalamus -> L4 and L6 -> L5 but
   has no corticothalamic feedback and no reticulo-thalamic inhibition, so it
   is not a loop. We add ``L6_E -> TCR``/``L6_E -> nRT`` (corticothalamic
   feedback) and ``nRT -> TCR``/``nRT -> nRT`` (the reciprocal reticular
   inhibition). The TCR <-> nRT reciprocal loop is the canonical spindle
   generator; its delays are tuned so the round trip resonates near 13 Hz.

2. **Sleep drives.** A ~1 Hz ``ac_generator`` injects current into the cortical
   pyramidal pools (and a matched slow modulation into thalamus) to cycle the
   network between depolarised UP states and silent DOWN states (the slow
   oscillation). A ~13 Hz ``ac_generator`` into nRT/TCR, made spindle-permissive
   only on the UP phase, produces waxing/waning spindles nested on the slow
   wave.

Honest scope: ``iaf_cond_exp`` has no intrinsic Ih / T-type Ca2+ currents, so
spindles here are an *imposed-drive + loop-resonance* hybrid rather than purely
intrinsic thalamic rhythmogenesis. Swapping the neuron model to NEST's
``ht_neuron`` (Hill-Tononi) would make both rhythms emergent; this is left as a
documented upgrade for the bio-plausible runs.

The thalamus is interpreted as the **MGB (medial geniculate body)**, the
auditory relay; its drive represents auditory-nerve / inferior-colliculus input.
"""

import math
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np


# ---------------------------------------------------------------------------
#  Configuration dataclasses (adapted from cc core/simulator.py)
# ---------------------------------------------------------------------------


@dataclass
class NetworkConfig:
    """Population sizes, subtype splits, and integration parameters."""

    N_thalamus_E: int = 5   # TCR / MGB relay cells
    N_thalamus_I: int = 1   # nRT reticular cells
    N_L4_E: int = 24
    N_L4_I: int = 6
    N_L23_E: int = 24
    N_L23_I: int = 6
    N_L5_E: int = 20
    N_L5_I: int = 5
    N_L6_E: int = 16
    N_L6_I: int = 4
    split_L23_E: int = 2
    split_L23_I: int = 3
    split_L4_E: int = 1
    split_L4_I: int = 1
    split_L5_E: int = 2
    split_L56_I: int = 3
    split_L6_E: int = 1
    tstop: float = 3000.0
    dt: float = 0.1
    v_init: float = -65.0

    @classmethod
    def from_dict(cls, config: dict) -> "NetworkConfig":
        kwargs: Dict = {}
        net = config.get("network", {})
        for layer, prefix in [("thalamus", "N_thalamus"), ("L4", "N_L4"),
                              ("L23", "N_L23"), ("L5", "N_L5"), ("L6", "N_L6")]:
            if layer in net:
                kwargs[f"{prefix}_E"] = net[layer].get("E", getattr(cls, f"{prefix}_E"))
                kwargs[f"{prefix}_I"] = net[layer].get("I", getattr(cls, f"{prefix}_I"))
        for key in ["L23_E", "L23_I", "L4_E", "L4_I", "L5_E", "L56_I", "L6_E"]:
            if key in config.get("splits", {}):
                kwargs[f"split_{key}"] = config["splits"][key]
        sim = config.get("simulation", {})
        for key in ["tstop", "dt", "v_init"]:
            if key in sim:
                kwargs[key] = sim[key]
        return cls(**kwargs)

    @classmethod
    def from_file(cls, path: Union[str, Path]) -> "NetworkConfig":
        path = Path(path)
        if not path.exists():
            print(f"Config file not found: {path}. Using defaults.")
            return cls()
        with open(path, "r") as f:
            if path.suffix in (".yaml", ".yml"):
                try:
                    import yaml
                    config = yaml.safe_load(f)
                except ImportError:
                    raise SystemExit(
                        "PyYAML is required to read YAML configs "
                        "(`pip install pyyaml`), or pass a .json config."
                    )
            elif path.suffix == ".json":
                config = json.load(f)
            else:
                raise ValueError(f"Unsupported config format: {path.suffix}")
        return cls.from_dict(config)

    @property
    def total_neurons(self) -> int:
        return (self.N_thalamus_E + self.N_thalamus_I + self.N_L4_E + self.N_L4_I
                + self.N_L23_E + self.N_L23_I + self.N_L5_E + self.N_L5_I
                + self.N_L6_E + self.N_L6_I)

    def summary(self) -> str:
        return "\n".join([
            "Auditory thalamo-cortical network:",
            f"  MGB/thalamus: {self.N_thalamus_E}E (TCR) / {self.N_thalamus_I}I (nRT)",
            f"  L4:           {self.N_L4_E}E / {self.N_L4_I}I",
            f"  L2/3:         {self.N_L23_E}E / {self.N_L23_I}I",
            f"  L5:           {self.N_L5_E}E / {self.N_L5_I}I",
            f"  L6:           {self.N_L6_E}E / {self.N_L6_I}I",
            f"  total:        {self.total_neurons} neurons, "
            f"tstop={self.tstop} ms, dt={self.dt} ms",
        ])


@dataclass
class SynapseParams:
    """Synaptic weights/delays and connection probability.

    Connection *topology* and weight regime follow cc; the synaptic
    **biophysics** (reversal potentials and decay time constants) are taken from
    the conductance-based thalamocortical sleep model of Mushtaq et al. 2024
    (PLOS ONE 19(6):e0306218), Table 3 and the "Synaptic currents" section.

    Note on units: Mushtaq's maximal conductances are HH area-based
    (mS/cm^2 x area), so the *absolute* values are not transferable to our point
    ``iaf_cond_exp`` cells; what transfers are the reversal potentials, the
    kinetics, and the relative conductance ratios (see ``SleepParams``).
    """

    # Charge-preserving calibration: adopting the slower Mushtaq AMPA tau
    # (2.0 -> 5.3 ms) raises the charge per excitatory event ~2.6x, which alone
    # drives the recurrent cortex into runaway. iaf_cond_exp cannot inherit the
    # HH area-based absolute conductances, so the intracortical excitatory
    # weight is scaled by ~2.0/5.3 to hold weight x tau (the per-event charge)
    # roughly constant and keep cortex in a plausible firing regime.
    exc_weight_mean: float = 0.00038  # uS (NEURON ExpSyn units; scaled to nS)
    exc_weight_std: float = 0.00034
    # AMPA decay tau. Mushtaq 2024 first-order scheme: alpha=1.1, beta=0.19 ms^-1
    # -> tau_decay = 1/beta ~= 5.3 ms (NMDA, tau ~150 ms, is not represented by
    # the single exp conductance of iaf_cond_exp and is omitted).
    exc_tau: float = 5.3
    exc_delay_mean: float = 3.0
    exc_delay_std: float = 1.0
    exc_e: float = 0.0               # AMPA/NMDA E_syn = 0 mV (Mushtaq 2024)
    min_delay: float = 0.2
    inh_weight: float = 0.0015
    # GABA_A decay tau. Mushtaq 2024: beta=0.166 ms^-1 -> tau_decay ~= 6.0 ms.
    inh_tau: float = 6.0
    inh_delay: float = 2.0
    inh_e: float = -70.0             # GABA_A E_syn = -70 mV (Mushtaq 2024)
    conn_prob: float = 0.1


@dataclass
class SleepParams:
    """Drives that produce the sleep slow-wave and spindles.

    The slow oscillation is an injected ~1 Hz current that swings cortex (and,
    weaker, thalamus) between UP and DOWN. Spindles are a ~13 Hz current into
    the reticular/relay cells that only entrains during the UP phase, reinforced
    by the tuned nRT<->TCR reciprocal loop.
    """

    # auditory-nerve / brainstem background drive onto MGB (Hz, per relay cell).
    # Kept very light: with conductance synapses the driving force at the relay's
    # operating point (~ -57 mV to E_ex = 0) is large, so even a few nS of tonic
    # background saturates the relay and erases the spindle rhythm. During sleep
    # the relay is rhythmic, not in tonic relay mode, so this is also biological.
    auditory_rate: float = 20.0       # Hz
    auditory_weight: float = 0.0006   # uS

    # slow oscillation (~1 Hz). With iaf_cond_exp here R ~= 100 MOhm, i.e.
    # ~0.1 mV per pA; V_th-E_L is 15 mV, so ~150 pA reaches threshold.
    slow_freq: float = 1.0            # Hz
    slow_amp_cortex: float = 120.0    # pA
    slow_offset_cortex: float = 50.0  # pA (DC bias toward threshold)
    # Thalamus operating point (rheobase ~= 150 pA): on the UP phase the DC
    # baseline (offset+amp ~= 150 pA) sits right at threshold so the 13 Hz
    # spindle peaks make it fire ~13-15 Hz; on the DOWN trough the baseline
    # drops to ~70 pA and even spindle peaks stay sub-threshold -> UP-gated.
    slow_amp_thalamus: float = 40.0   # pA
    slow_offset_thalamus: float = 110.0  # pA

    # spindle (~13 Hz): peaks push the (UP-phase, near-threshold) relay/reticular
    # cells over threshold once per cycle -> 13 Hz-locked firing nested on UP.
    spindle_freq: float = 13.0        # Hz
    spindle_amp: float = 60.0         # pA, into nRT + TCR
    spindle_phase_deg: float = 0.0

    # Reticulo-thalamic loop tuning (the spindle resonator). The connection
    # topology and the *relative* weighting follow Mushtaq et al. 2024 Table 3
    # (intra-thalamic block): TC->RE (AMPA, .1 uS), RE->TC (GABA_A .05 + GABA_B
    # .02 uS), RE->RE (GABA_A .05 uS). GABA_B is lumped into the single GABA
    # conductance here, so RE->TC is the strongest inhibition (~.07) and
    # RE->RE matches GABA_A (.05). Absolute uS are not transferable (HH area-
    # based -> point cell), so these are calibrated to our operating point while
    # keeping the resonator at ~13 Hz.
    rt_inh_weight: float = 0.014      # uS, nRT -> TCR  (Mushtaq RE->TC, strongest)
    rt_inh_delay: float = 6.0         # ms (half of ~1/13 s round trip)
    tr_exc_weight: float = 0.006      # uS, TCR -> nRT  (Mushtaq TC->RE, AMPA)
    tr_exc_delay: float = 4.0         # ms
    nrt_self_weight: float = 0.006    # uS, nRT -> nRT  (Mushtaq RE->RE, GABA_A)
    # Corticothalamic feedback (L6 -> thalamus), closes the loop. Mushtaq 2024
    # Table 3 (cortico-thalamic block) has PY->TC (AMPA .003) twice PY->RE
    # (AMPA .0015), so the feedback onto the relay is twice that onto nRT.
    ct_fb_weight: float = 0.0008      # uS, L6_E -> TCR  (Mushtaq PY->TC)
    ct_fb_weight_nrt: float = 0.0004  # uS, L6_E -> nRT  (Mushtaq PY->RE, = TC/2)
    # Slow GABA_B component of RE->TC (Mushtaq 2024 Table 3, .02 uS). Only used
    # by the ht_neuron (HH) model; its long time constant paces the spindle
    # envelope and de-inactivates I_T for the rebound burst.
    rt_inh_gaba_b: float = 0.02       # uS, nRT -> TCR (GABA_B; ht_neuron only)
    # HH (ht_neuron) emergent mode: drop the imposed 13 Hz spindle AC; spindles
    # then emerge from the RE<->TC loop (I_T rebound + GABA_B). The light 1 Hz AC
    # is kept as a scaffold that gates the emergent spindles to UP states.
    emergent_spindles: bool = False

    # --- external modulatory spindle trigger (HH mode; Fernandez & Luthi 2020)
    # Design rule: the modulatory signal decides *when* a spindle is initiated
    # and *whether* it is permitted, but never its frequency -- the intra-spindle
    # rhythm stays emergent from the RE<->TC loop. See
    # docs/spindle_review_mapping.md sect. 5.
    #
    # (A) Phasic trigger: a brief depolarising "kick" onto RE/nRT once per slow
    # oscillation, standing in for the layer-6 corticothalamic burst that
    # initiates spindles (review sect. V.B.2).
    spindle_trigger: bool = False
    trigger_amp: float = 30.0        # pA into nRT
    trigger_dur: float = 30.0        # ms pulse width
    trigger_phase_ms: float = 250.0  # ms after cycle start (SO UP onset)
    # Minimum interval between trigger pulses, snapped to a whole number of SO
    # cycles so spindles stay UP-state-locked. 0 = fire every cycle (default).
    #
    # NOTE (measured): spacing the trigger does NOT by itself produce the
    # review's 5-10 s spindle refractory -- with the trigger at 7 s the model
    # still yields ~1 spindle/s, because the 1 Hz SO drive plus L6
    # corticothalamic feedback re-activate the thalamus on every UP state
    # regardless of the trigger. See docs/spindle_review_mapping.md sect. 5.4.
    trigger_refractory_ms: float = 0.0


@dataclass
class HHParams:
    """Hodgkin-Huxley (``ht_neuron``, Hill & Tononi 2005) intrinsic parameters.

    ``ht_neuron`` carries the intrinsic currents needed for *emergent*
    thalamocortical sleep rhythms -- no imposed spindle oscillator:

    - **I_T** (low-threshold T-type Ca2+): de-inactivates when the cell is
      hyperpolarised, giving the post-inhibitory **rebound burst** that, round
      the RE<->TC loop, is the spindle.
    - **I_h** (hyperpolarisation-activated cation): the pacemaker/refractory
      current that makes spindles **wax and wane** and terminate.
    - **I_NaP / I_KNa** (persistent Na+ / Na-dependent K+): the cortical
      depolarising drive and spike-frequency adaptation behind UP/DOWN states.

    The **sleep vs wake regime is set by the K-leak conductance g_KL** (Hill &
    Tononi's neuromodulation knob): a larger g_KL hyperpolarises the cell toward
    the burst-permissive range. Thalamic cells get the larger g_KL so I_T
    de-inactivates and spindles can emerge.
    """

    # K-leak (neuromodulatory sleep depth). Larger -> more hyperpolarised.
    g_KL_cortex: float = 1.0
    g_KL_thalamus: float = 1.0     # burst-permissive but still able to fire
    g_NaL: float = 0.2
    # intrinsic current peak conductances (nS); ht_neuron defaults are sensible,
    # thalamic I_T / I_h are raised so the RE<->TC rebound loop actually bursts.
    g_peak_T_thalamus: float = 3.0
    g_peak_h_thalamus: float = 2.0
    # SK2 surrogate in the reticular (RE/nRT) cells. Fernandez & Luthi 2020
    # (sect. V.A.1): Ca2+ entering TRN dendrites through Ca_v3.3 activates SK2,
    # producing the burst after-hyperpolarisation that keeps TRN bursts short
    # and repeatable -- which shapes spindle amplitude and regularity.
    # ht_neuron has no Ca-gated K+ current, but its Na-dependent K+ current
    # (I_KNa) is likewise burst-activated and produces an AHP, so it stands in
    # for SK2 functionally (see docs/spindle_review_mapping.md sect. 4; the
    # faithful version is a NESTML ht_neuron_sk with a fast Ca2+ pool).
    g_peak_KNa_RE: float = 0.5
    # TC Ca2+ handling, which drives the I_h after-depolarisation behind spindle
    # TERMINATION and the 5-10 s refractory period (Fernandez & Luthi sect.
    # V.D.1: Ca2+ entry via Ca_v3.1 -> cAMP -> HCN, stabilising I_h open state).
    # ht_neuron models exactly this; beta_Ca sets Ca2+ influx per T-current and
    # tau_Ca its decay (default 10 s -- the right refractory timescale).
    beta_Ca_TC: float = 0.001
    tau_Ca_TC: float = 10000.0
    g_peak_NaP_cortex: float = 1.0
    g_peak_KNa_cortex: float = 1.0
    tau_m: float = 16.0
    theta_eq: float = -51.0
    # ht_neuron drive amplitudes (pA). ht_neuron's rheobase is ~35 pA (much
    # lower than the iaf 150 pA), so the imposed-iaf amplitudes would saturate
    # it. Cortex: a light 1 Hz scaffold on top of the emergent I_NaP/I_KNa
    # bistability. Thalamus: held just below threshold and lifted toward it on
    # the UP phase, so the RE<->TC loop rings (emergent spindle) UP-locked.
    cortex_slow_amp: float = 32.0
    cortex_slow_offset: float = 4.0
    # TC relay needs a stronger depolarising baseline than RE so it fires in the
    # gaps between RE inhibitory volleys (the spindle); RE only needs enough to
    # respond to corticothalamic feedback.
    tc_slow_amp: float = 12.0
    tc_slow_offset: float = 34.0
    thal_slow_amp: float = 12.0    # RE (reticular)
    thal_slow_offset: float = 14.0
    # ht_neuron synapses scale g_peak by the connection weight; the intra-
    # thalamic and corticothalamic loop weights need boosting (vs the iaf uS
    # regime) for the RE<->TC loop to actually engage and ring.
    loop_weight_gain: float = 3.0
    aud_weight_gain: float = 8.0
    # (B) Slow neuromodulatory "permission" signal (Fernandez & Luthi 2020,
    # sect. IV.E): sigma power / spindle occurrence is clustered on a ~0.02 Hz
    # (~50 s) infraslow rhythm that alternates NREMS between spindle-rich
    # CONTINUITY and spindle-poor FRAGILITY periods. Modelled as a slow swing of
    # thalamic excitability (the current-injection analogue of Hill & Tononi's
    # g_KL neuromodulation knob): on the permissive phase a triggered spindle
    # succeeds, on the opposite phase it does not -> clustering + refractoriness.
    # Note the swing must be LARGE: hyperpolarising the thalamus lowers TC but
    # also makes RE *more* burst-prone (I_T de-inactivates), so the two effects
    # partly cancel and only a big excursion modulates sigma power. Calibrated
    # by sweeping the sigma-envelope modulation ratio (see the runner's
    # "sigma infraslow mod" metric): ~40 pA gives a clear peak, ~20 pA none.
    infraslow_freq: float = 0.02     # Hz (~50 s period)
    infraslow_amp: float = 40.0      # pA swing of the thalamic baseline
    # Calibrated TOTAL intra-thalamic conductance onto each post-synaptic cell
    # (fan-in normalised; single-cell tests put the spindle-generating RE->TC
    # inhibition around these values). GABA_A + slow GABA_B for RE->TC.
    tot_re_tc_gaba_a: float = 45.0
    tot_re_tc_gaba_b: float = 22.0
    tot_tc_re_ampa: float = 30.0
    tot_re_re_gaba_a: float = 18.0


@dataclass
class SimulationConfig:
    num_threads: int = 1
    verbose: bool = False
    record_traces: bool = True
    seed: Optional[int] = None
    # "iaf_cond_exp" (default, point LIF) or "ht_neuron" (Hodgkin-Huxley,
    # Hill-Tononi, with intrinsic I_T/I_h/I_NaP/I_KNa -> emergent spindles).
    neuron_model: str = "iaf_cond_exp"


# ---------------------------------------------------------------------------
#  Simulator
# ---------------------------------------------------------------------------

# iaf_cond_exp gives conductance-based exponential synapses with separate
# E_ex / E_in reversal potentials, matching cc's choice and the NEURON ExpSyn.
_PREFERRED_NEURON_MODELS = ["iaf_cond_exp", "aeif_cond_exp"]


class AuditoryThalamoCorticalSleep:
    """Builds and runs the closed-loop auditory column under sleep drives."""

    WEIGHT_SCALE = 1000.0  # NEURON ExpSyn (uS) -> NEST conductance (nS)

    # Reference excitatory layer size (the local/cc default L23_E). Synaptic
    # weights are scaled by REF_E / N_L23_E so that, with fixed pairwise_bernoulli
    # probability, the expected per-neuron input current (and the E/I balance)
    # stays scale-invariant as the network grows for the bio-plausible MN5 run.
    # Without this, in-degree ~ p*N would scale with population size and a 10x
    # larger column would drive the cortex into saturation.
    REF_E = 24

    _nest_warning_shown = False

    def __init__(self, network_config: Optional[NetworkConfig] = None,
                 syn: Optional[SynapseParams] = None,
                 sleep: Optional[SleepParams] = None,
                 sim_config: Optional[SimulationConfig] = None):
        self.cfg = network_config or NetworkConfig()
        self.syn = syn or SynapseParams()
        self.sleep = sleep or SleepParams()
        self.sim = sim_config or SimulationConfig()

        # balanced-network weight scaling (1.0 for the local default sizes)
        self._w_scale = self.REF_E / max(1, self.cfg.N_L23_E)

        self._nest = None
        self._nest_available = False
        self._neuron_model: Optional[str] = None
        self._is_cond_model = False
        self._is_ht = False
        self.hh = HHParams()
        self._setup_nest()

    # -- setup -------------------------------------------------------------

    def _setup_nest(self):
        try:
            import nest
            self._nest = nest
            self._nest_available = True
        except ImportError:
            if not AuditoryThalamoCorticalSleep._nest_warning_shown:
                print("Note: NEST not available. The module imports, but run() "
                      "requires NEST 3.x.")
                AuditoryThalamoCorticalSleep._nest_warning_shown = True
            self._nest_available = False

    def _init_kernel(self):
        nest = self._nest
        nest.ResetKernel()
        resolution = max(0.025, self.cfg.dt)
        kernel = {
            "resolution": resolution,
            "local_num_threads": max(1, self.sim.num_threads),
            "print_time": self.sim.verbose,
        }
        if self.sim.seed is not None:
            kernel["rng_seed"] = self.sim.seed
        nest.set(**kernel)
        self._neuron_model = self._select_neuron_model()
        self._is_cond_model = "cond" in self._neuron_model
        self._is_ht = self._neuron_model == "ht_neuron"
        if self._is_ht:
            self._ht_receptors = self._nest.GetDefaults("ht_neuron")["receptor_types"]

    def _select_neuron_model(self) -> str:
        nest = self._nest
        # honour an explicit request (e.g. "ht_neuron" for the HH model)
        requested = getattr(self.sim, "neuron_model", None)
        if requested:
            if requested in nest.node_models:
                return requested
            print(f"Requested neuron model {requested!r} unavailable; "
                  f"falling back to the preferred list.")
        for model in _PREFERRED_NEURON_MODELS:
            if model in nest.node_models:
                return model
        raise RuntimeError(f"No suitable neuron model; tried {_PREFERRED_NEURON_MODELS}")

    def _neuron_params(self, role: str = "cortex") -> dict:
        """Per-role neuron parameters. ``role`` in {"cortex", "TC", "RE"}.

        For ``ht_neuron`` the role selects which intrinsic currents are active,
        mirroring Mushtaq 2024 Tables 1-2: cortex carries persistent-Na (I_NaP)
        and adaptation (I_KNa) but no I_T/I_h; thalamic relay (TC) carries I_T
        and I_h; reticular (RE) carries I_T only. Thalamic cells also get a
        larger K-leak (g_KL) so they rest hyperpolarised and I_T de-inactivates.
        """
        if self._is_ht:
            return self._ht_neuron_params(role)
        p = {
            "V_m": self.cfg.v_init,
            "tau_syn_ex": max(0.1, self.syn.exc_tau),
            "tau_syn_in": max(0.1, self.syn.inh_tau),
            "E_L": -70.0,
            "V_th": -55.0,
            "V_reset": -60.0,
            "C_m": 100.0,
            "t_ref": 2.0,
        }
        if self._is_cond_model:
            p.update({"E_ex": self.syn.exc_e, "E_in": self.syn.inh_e, "g_L": 10.0})
        return p

    def _ht_neuron_params(self, role: str) -> dict:
        hh = self.hh
        thalamic = role in ("TC", "RE")
        p = {
            "V_m": -70.0,
            "g_NaL": hh.g_NaL,
            "g_KL": hh.g_KL_thalamus if thalamic else hh.g_KL_cortex,
            "tau_m": hh.tau_m,
            "theta_eq": hh.theta_eq,
            "theta": hh.theta_eq,
            # receptor kinetics from the (Mushtaq-aligned) synapse params
            "tau_decay_AMPA": max(0.5, self.syn.exc_tau),
            "tau_decay_GABA_A": max(0.5, self.syn.inh_tau),
            "instant_unblock_NMDA": True,
        }
        if thalamic:
            # I_T (both) and I_h (relay only) on; no cortical NaP.
            # RE additionally gets the SK2-surrogate AHP current (I_KNa).
            p.update({
                "g_peak_NaP": 0.0,
                "g_peak_KNa": hh.g_peak_KNa_RE if role == "RE" else 0.0,
                "g_peak_T": hh.g_peak_T_thalamus,
                "g_peak_h": hh.g_peak_h_thalamus if role == "TC" else 0.0,
            })
            if role == "TC":
                # Ca2+ -> I_h pathway: the spindle-termination / refractory
                # mechanism (see HHParams.beta_Ca_TC).
                p.update({"beta_Ca": hh.beta_Ca_TC, "tau_Ca": hh.tau_Ca_TC})
        else:
            # cortical persistent-Na drive + Na-dependent-K adaptation; no T/h
            p.update({
                "g_peak_NaP": hh.g_peak_NaP_cortex,
                "g_peak_KNa": hh.g_peak_KNa_cortex,
                "g_peak_T": 0.0,
                "g_peak_h": 0.0,
            })
        return p

    # -- population creation ----------------------------------------------

    def _create_split(self, n_total: int, n_splits: int, params: dict) -> list:
        nest = self._nest
        if n_total <= 0:
            return []
        nodes = nest.Create(self._neuron_model, n_total, params=params)
        if n_splits <= 1:
            return [nodes]
        step = max(1, n_total // n_splits)
        out = []
        for i in range(n_splits):
            start = i * step
            end = n_total if i == n_splits - 1 else min(n_total, (i + 1) * step)
            if start < n_total:
                out.append(nodes[start:end])
        return out

    def build_network(self) -> Dict[str, list]:
        cfg = self.cfg
        np_ = self._neuron_params("cortex")
        np_tc = self._neuron_params("TC")   # relay: I_T + I_h
        np_re = self._neuron_params("RE")   # reticular: I_T only
        n56i = cfg.N_L5_I + cfg.N_L6_I
        return {
            "thalamus_E": [self._nest.Create(self._neuron_model, cfg.N_thalamus_E, params=np_tc)],
            "thalamus_I": [self._nest.Create(self._neuron_model, cfg.N_thalamus_I, params=np_re)],
            "L4_E":       self._create_split(cfg.N_L4_E, cfg.split_L4_E, np_),
            "L4_I":       self._create_split(cfg.N_L4_I, cfg.split_L4_I, np_),
            "L23_E_RS":   self._create_split(cfg.N_L23_E, cfg.split_L23_E, np_),
            "L23_E_FRB":  self._create_split(cfg.N_L23_E, cfg.split_L23_E, np_),
            "L23_I_Bask": self._create_split(cfg.N_L23_I, cfg.split_L23_I, np_),
            "L23_I_LTS":  self._create_split(cfg.N_L23_I, cfg.split_L23_I, np_),
            "L23_I_Axax": self._create_split(cfg.N_L23_I, cfg.split_L23_I, np_),
            "L5_E_RS":    self._create_split(cfg.N_L5_E, cfg.split_L5_E, np_),
            "L5_E_IB":    self._create_split(cfg.N_L5_E, cfg.split_L5_E, np_),
            "L56_I_Bask": self._create_split(n56i, cfg.split_L56_I, np_),
            "L56_I_LTS":  self._create_split(n56i, cfg.split_L56_I, np_),
            "L56_I_Axax": self._create_split(n56i, cfg.split_L56_I, np_),
            "L6_E":       self._create_split(cfg.N_L6_E, cfg.split_L6_E, np_),
        }

    # -- connection helpers (adapted from cc) -----------------------------

    @staticmethod
    def flatten(population):
        if isinstance(population, list):
            if not population:
                return None
            result = population[0]
            for nc in population[1:]:
                result = result + nc
            return result
        return population

    @staticmethod
    def _merge(*node_collections):
        valid = [nc for nc in node_collections if nc is not None and len(nc) > 0]
        if not valid:
            return None
        result = valid[0]
        for nc in valid[1:]:
            result = result + nc
        return result

    def _min_delay(self) -> float:
        return max(max(0.025, self.cfg.dt), self.syn.min_delay)

    def _ht_gain(self, w):
        """Boost a loop weight for ht_neuron (g_peak-scaled) synapses; identity
        for iaf_cond_exp so its validated tuning is untouched."""
        return w * self.hh.loop_weight_gain if self._is_ht else w

    def _ht_thal_w(self, iaf_weight, target_total, n_presyn):
        """Weight for an all-to-all (prob=1.0) intra-thalamic loop synapse.

        iaf_cond_exp: return its validated uS weight unchanged. ht_neuron: pick
        a per-synapse weight so the TOTAL conductance onto each post-synaptic
        cell (summed over all ``n_presyn`` inputs) equals a calibrated
        ``target_total`` -- fan-in normalisation, so the RE<->TC loop behaves the
        same whether the thalamus has 4 or 40 reticular cells (without it, 40 RE
        cells crush every TC and no rebound spindle can form)."""
        if not self._is_ht:
            return iaf_weight
        return target_total / (max(1, n_presyn) * self.WEIGHT_SCALE * self._w_scale)

    def _exc_syn_spec(self, weight_nS, delay_ms) -> dict:
        if self._is_ht:
            # ht_neuron routes by receptor_type -> AMPA. Excitatory weights are
            # already positive (no abs(): weight_nS may be a NEST Parameter).
            return {"synapse_model": "static_synapse", "weight": weight_nS,
                    "delay": delay_ms, "receptor_type": self._ht_receptors["AMPA"]}
        # iaf_cond_exp / aeif_cond_exp route by weight SIGN (positive -> g_ex /
        # E_ex), with no receptor_type. Keep weights positive for excitation.
        return {"synapse_model": "static_synapse", "weight": weight_nS, "delay": delay_ms}

    def _inh_syn_spec(self, weight_nS, delay_ms, receptor="GABA_A") -> dict:
        w = abs(weight_nS)
        if self._is_ht:
            # ht_neuron: positive weight into the GABA_A (or GABA_B) receptor.
            return {"synapse_model": "static_synapse", "weight": w,
                    "delay": delay_ms, "receptor_type": self._ht_receptors[receptor]}
        # Inhibition = NEGATIVE weight (routes to g_in / E_in). Critically, NOT
        # a positive weight + receptor_type: cond models reject receptor_type,
        # which would silently turn inhibition into excitation (disinhibition).
        return {"synapse_model": "static_synapse", "weight": -w, "delay": delay_ms}

    def _safe_connect(self, source, target, conn_spec, syn_spec):
        self._nest.Connect(source, target, conn_spec=conn_spec, syn_spec=syn_spec)

    def connect_exc(self, source_pop, target_pop, weight=None, delay=None, prob=None):
        source = self.flatten(source_pop)
        target = self.flatten(target_pop)
        if source is None or target is None or len(source) == 0 or len(target) == 0:
            return
        p = self.syn
        w_mean = (weight if weight is not None else p.exc_weight_mean) * self.WEIGHT_SCALE * self._w_scale
        w_std = p.exc_weight_std * self.WEIGHT_SCALE * self._w_scale
        min_d = self._min_delay()
        d_mean = max(min_d, delay if delay is not None else p.exc_delay_mean)

        if weight is None and w_std > 0 and w_mean > 0:
            sigma_ln = math.sqrt(math.log(1 + (w_std / w_mean) ** 2))
            mu_ln = math.log(w_mean) - sigma_ln ** 2 / 2
            weight_param = self._nest.random.lognormal(mean=mu_ln, std=sigma_ln)
            d_lo = max(min_d, d_mean - 2 * p.exc_delay_std)
            d_hi = max(d_lo + 0.1, d_mean + 2 * p.exc_delay_std)
            delay_param = self._nest.random.uniform(min=d_lo, max=d_hi)
        else:
            weight_param = w_mean
            delay_param = d_mean

        syn_spec = self._exc_syn_spec(weight_param, delay_param)
        conn_spec = {"rule": "pairwise_bernoulli",
                     "p": prob if prob is not None else p.conn_prob}
        self._safe_connect(source, target, conn_spec, syn_spec)

    def connect_inh(self, source_pop, target_pop, weight=None, delay=None,
                    prob=None, gaba_b=None):
        source = self.flatten(source_pop)
        target = self.flatten(target_pop)
        if source is None or target is None or len(source) == 0 or len(target) == 0:
            return
        p = self.syn
        w = (weight if weight is not None else p.inh_weight) * self.WEIGHT_SCALE * self._w_scale
        min_d = self._min_delay()
        d = max(min_d, delay if delay is not None else p.inh_delay)
        conn_spec = {"rule": "pairwise_bernoulli",
                     "p": prob if prob is not None else p.conn_prob}
        self._safe_connect(source, target, conn_spec, self._inh_syn_spec(w, d))
        # ht_neuron only: an additional slow GABA_B component (e.g. RE->TC).
        # GABA_B's long time constant paces the waxing/waning of the spindle.
        if self._is_ht and gaba_b:
            wb = gaba_b * self.WEIGHT_SCALE * self._w_scale
            self._safe_connect(source, target, conn_spec,
                               self._inh_syn_spec(wb, d, receptor="GABA_B"))

    # -- network wiring ----------------------------------------------------

    def connect_network(self, pops: Dict[str, list]) -> None:
        """Feed-forward + recurrent cortical wiring (cc), plus the closed
        thalamo-cortical loop additions."""
        # --- feed-forward / intracortical excitation (cc topology) ---
        # MGB -> L4 is the dense thalamocortical input; with few MGB cells we
        # raise prob/weight so L4 is actually driven (and spindles propagate).
        self.connect_exc(pops["thalamus_E"], pops["L4_E"],
                         weight=self._ht_gain(0.004),
                         delay=self.syn.exc_delay_mean, prob=0.6)
        # Thalamocortical feedforward inhibition: TC -> IN (Mushtaq 2024 Table 3,
        # TC->IN AMPA, same conductance as TC->PY but a smaller connecting radius
        # (5 vs 21) -> more focal). This was missing in cc; it lets thalamic input
        # recruit cortical interneurons, sharpening the UP-state response.
        self.connect_exc(pops["thalamus_E"], pops["L4_I"],
                         weight=0.004, delay=self.syn.exc_delay_mean, prob=0.3)
        self.connect_exc(pops["thalamus_E"], pops["thalamus_I"],  # TCR -> nRT
                         weight=self._ht_thal_w(self.sleep.tr_exc_weight,
                                                self.hh.tot_tc_re_ampa,
                                                self.cfg.N_thalamus_E),
                         delay=self.sleep.tr_exc_delay, prob=1.0)
        self.connect_exc(pops["L4_E"], pops["L4_E"])
        self.connect_exc(pops["L4_E"], pops["L23_E_RS"])
        self.connect_exc(pops["L4_E"], pops["L23_E_FRB"])
        self.connect_exc(pops["L4_E"], pops["L5_E_RS"])
        self.connect_exc(pops["L4_E"], pops["L5_E_IB"])
        self.connect_exc(pops["L4_E"], pops["L4_I"])

        self.connect_exc(pops["L23_E_RS"], pops["L23_E_RS"])
        self.connect_exc(pops["L23_E_FRB"], pops["L23_E_FRB"])
        self.connect_exc(pops["L23_E_RS"], pops["L23_E_FRB"])
        self.connect_exc(pops["L23_E_FRB"], pops["L23_E_RS"])
        self.connect_exc(pops["L23_E_RS"], pops["L5_E_RS"])
        self.connect_exc(pops["L23_E_RS"], pops["L5_E_IB"])
        self.connect_exc(pops["L23_E_RS"], pops["L6_E"])
        self.connect_exc(pops["L23_E_FRB"], pops["L5_E_RS"])
        self.connect_exc(pops["L23_E_FRB"], pops["L5_E_IB"])
        self.connect_exc(pops["L23_E_FRB"], pops["L6_E"])

        self.connect_exc(pops["L5_E_RS"], pops["L5_E_RS"])
        self.connect_exc(pops["L5_E_IB"], pops["L5_E_IB"])
        self.connect_exc(pops["L5_E_RS"], pops["L5_E_IB"])
        self.connect_exc(pops["L5_E_IB"], pops["L5_E_RS"])
        self.connect_exc(pops["L5_E_RS"], pops["L6_E"])
        self.connect_exc(pops["L5_E_IB"], pops["L6_E"])
        self.connect_exc(pops["L6_E"], pops["L5_E_RS"])
        self.connect_exc(pops["L6_E"], pops["L5_E_IB"])

        # --- intracortical inhibition (cc topology) ---
        self.connect_inh(pops["L4_I"], pops["L4_E"])
        self.connect_inh(pops["L23_I_LTS"], pops["L23_E_RS"])
        self.connect_inh(pops["L23_I_LTS"], pops["L23_E_FRB"])
        self.connect_inh(pops["L23_I_Bask"], pops["L23_E_RS"])
        self.connect_inh(pops["L23_I_Bask"], pops["L23_E_FRB"])
        self.connect_inh(pops["L56_I_LTS"], pops["L5_E_RS"])
        self.connect_inh(pops["L56_I_LTS"], pops["L5_E_IB"])
        self.connect_inh(pops["L56_I_Bask"], pops["L5_E_RS"])
        self.connect_inh(pops["L56_I_Bask"], pops["L6_E"])

        # --- CLOSING THE LOOP (additions over cc) ---
        s = self.sleep
        # corticothalamic feedback: L6 -> MGB (TCR) and L6 -> nRT. Mushtaq 2024
        # Table 3: PY->TC (.003) is twice PY->RE (.0015), so the relay gets the
        # stronger feedback (ct_fb_weight) and nRT the weaker (ct_fb_weight_nrt).
        self.connect_exc(pops["L6_E"], pops["thalamus_E"],
                         weight=self._ht_gain(s.ct_fb_weight), delay=s.tr_exc_delay)
        self.connect_exc(pops["L6_E"], pops["thalamus_I"],
                         weight=self._ht_gain(s.ct_fb_weight_nrt), delay=s.tr_exc_delay)
        # reciprocal reticulo-thalamic inhibition: nRT -> TCR and nRT -> nRT.
        # The TCR<->nRT loop with these delays is the spindle resonator. Under
        # ht_neuron, RE->TC also gets the slow GABA_B component (Mushtaq 2024
        # Table 3, .02) whose long tau sets the spindle waxing/waning and,
        # crucially, hyperpolarises TC enough to de-inactivate I_T for the
        # rebound burst -- so the spindle is emergent, not imposed.
        n_re = self.cfg.N_thalamus_I
        self.connect_inh(pops["thalamus_I"], pops["thalamus_E"],
                         weight=self._ht_thal_w(s.rt_inh_weight,
                                                self.hh.tot_re_tc_gaba_a, n_re),
                         delay=s.rt_inh_delay, prob=1.0,
                         gaba_b=self._ht_thal_w(s.rt_inh_gaba_b,
                                                self.hh.tot_re_tc_gaba_b, n_re))
        self.connect_inh(pops["thalamus_I"], pops["thalamus_I"],
                         weight=self._ht_thal_w(s.nrt_self_weight,
                                                self.hh.tot_re_re_gaba_a, n_re),
                         delay=s.rt_inh_delay, prob=1.0)

    # -- drives ------------------------------------------------------------

    def attach_auditory_input(self, pops):
        """Auditory-nerve / brainstem background onto MGB relay cells."""
        nest = self._nest
        mgb = self.flatten(pops["thalamus_E"])
        if mgb is None or len(mgb) == 0:
            return
        gen = nest.Create("poisson_generator", len(mgb),
                          params={"rate": self.sleep.auditory_rate})
        aud_w = self.sleep.auditory_weight
        if self._is_ht:
            aud_w *= self.hh.aud_weight_gain   # ht g_peak-scaled synapses
        w = aud_w * self.WEIGHT_SCALE
        self._safe_connect(gen, mgb, {"rule": "one_to_one"},
                           self._exc_syn_spec(w, self._min_delay()))

    def _attach_spindle_trigger(self, re_pop):
        """Phasic corticothalamic-like kick onto RE, once per slow-oscillation
        cycle, standing in for the layer-6 burst that *initiates* a spindle.

        This sets only the initiation time: the spindle's frequency, duration
        and waxing/waning still emerge from the RE<->TC loop, and whether the
        kick actually produces a spindle depends on the slow (~0.02 Hz)
        permission signal. See docs/spindle_review_mapping.md sect. 5.
        """
        nest = self._nest
        s = self.sleep
        if re_pop is None or len(re_pop) == 0 or not s.trigger_amp:
            return
        so_period = 1000.0 / max(1e-6, s.slow_freq)   # ms per SO cycle
        # Fire only every n-th SO cycle so the inter-spindle interval respects
        # the 5-10 s refractory, while staying locked to an SO UP state.
        n_cycles = max(1, int(round(max(so_period, s.trigger_refractory_ms)
                                    / so_period)))
        step = n_cycles * so_period
        dur = max(self.cfg.dt, s.trigger_dur)
        times, amps = [], []
        t = max(self.cfg.dt, s.trigger_phase_ms)
        while t + dur < self.cfg.tstop:
            times += [t, t + dur]
            amps += [s.trigger_amp, 0.0]
            t += step
        if not times:
            return
        gen = nest.Create("step_current_generator", params={
            "amplitude_times": times, "amplitude_values": amps})
        nest.Connect(gen, re_pop)

    def attach_sleep_drives(self, pops):
        """Slow-wave (1 Hz) and spindle (13 Hz) current injections.

        - Cortex gets a 1 Hz sinusoid (+DC) so it cycles UP/DOWN -> slow wave.
        - Thalamus gets a weaker, *anti-phase-offset* 1 Hz term that lifts it
          to spindle-permissive only during the cortical UP phase, plus a 13 Hz
          sinusoid; the nRT<->TCR loop turns that into waxing/waning spindles.
        """
        nest = self._nest
        s = self.sleep

        cortex_E = self._merge(
            self.flatten(pops["L4_E"]),
            self.flatten(pops["L23_E_RS"]), self.flatten(pops["L23_E_FRB"]),
            self.flatten(pops["L5_E_RS"]), self.flatten(pops["L5_E_IB"]),
            self.flatten(pops["L6_E"]))
        thal = self._merge(self.flatten(pops["thalamus_E"]),
                           self.flatten(pops["thalamus_I"]))
        tc = self.flatten(pops["thalamus_E"])
        re = self.flatten(pops["thalamus_I"])

        emergent = self._is_ht and s.emergent_spindles

        # 1 Hz slow oscillation into cortex (UP/DOWN states). HH uses lighter
        # amplitudes (its rheobase is far lower than iaf's).
        cx_amp = self.hh.cortex_slow_amp if self._is_ht else s.slow_amp_cortex
        cx_off = self.hh.cortex_slow_offset if self._is_ht else s.slow_offset_cortex
        if cortex_E is not None and len(cortex_E) > 0:
            slow_cortex = nest.Create("ac_generator", params={
                "amplitude": cx_amp,
                "offset": cx_off,
                "frequency": s.slow_freq,
                "phase": 0.0,
            })
            nest.Connect(slow_cortex, cortex_E)

        if thal is not None and len(thal) > 0:
            if emergent:
                # HH emergent mode: no imposed 13 Hz oscillator. Light 1 Hz terms
                # hold TC/RE near threshold and lift them on the UP phase; the
                # cortical UP state (via L6->RE feedback) then kicks the RE<->TC
                # loop into self-generated, waxing/waning spindle bursts locked
                # to UP states. TC and RE are driven separately (TC needs the
                # stronger baseline to fire between RE inhibitory volleys); the
                # spindle frequency is set by the loop, not by any drive.
                if tc is not None and len(tc) > 0:
                    tc_drive = nest.Create("ac_generator", params={
                        "amplitude": self.hh.tc_slow_amp,
                        "offset": self.hh.tc_slow_offset,
                        "frequency": s.slow_freq, "phase": 0.0})
                    nest.Connect(tc_drive, tc)
                if re is not None and len(re) > 0:
                    re_drive = nest.Create("ac_generator", params={
                        "amplitude": self.hh.thal_slow_amp,
                        "offset": self.hh.thal_slow_offset,
                        "frequency": s.slow_freq, "phase": 0.0})
                    nest.Connect(re_drive, re)
                # (B) slow neuromodulatory permission signal: ~0.02 Hz swing of
                # thalamic excitability -> spindle-rich CONTINUITY vs
                # spindle-poor FRAGILITY periods (Fernandez & Luthi sect. IV.E).
                if self.hh.infraslow_amp:
                    infra = nest.Create("ac_generator", params={
                        "amplitude": self.hh.infraslow_amp,
                        "offset": 0.0,
                        "frequency": self.hh.infraslow_freq,
                        "phase": 0.0})
                    nest.Connect(infra, thal)
                # (A) phasic corticothalamic-like trigger onto RE
                if s.spindle_trigger:
                    self._attach_spindle_trigger(re)
            else:
                # imposed-drive mode (iaf_cond_exp): 1 Hz gating + 13 Hz spindle
                slow_thal = nest.Create("ac_generator", params={
                    "amplitude": s.slow_amp_thalamus,
                    "offset": s.slow_offset_thalamus,
                    "frequency": s.slow_freq,
                    "phase": 0.0,   # in-phase: thalamus depolarised during UP
                })
                nest.Connect(slow_thal, thal)

                spindle = nest.Create("ac_generator", params={
                    "amplitude": s.spindle_amp,
                    "offset": 0.0,
                    "frequency": s.spindle_freq,
                    "phase": s.spindle_phase_deg,
                })
                nest.Connect(spindle, thal)

    # -- recording ---------------------------------------------------------

    def _layer_nodes(self, pops):
        return {
            "MGB":  self.flatten(pops["thalamus_E"]),
            "nRT":  self.flatten(pops["thalamus_I"]),
            "L4":   self.flatten(pops["L4_E"]),
            "L23":  self._merge(self.flatten(pops["L23_E_RS"]), self.flatten(pops["L23_E_FRB"])),
            "L5":   self._merge(self.flatten(pops["L5_E_RS"]), self.flatten(pops["L5_E_IB"])),
            "L6":   self.flatten(pops["L6_E"]),
        }

    def _setup_recorders(self, pops):
        nest = self._nest
        layer_nodes = self._layer_nodes(pops)
        recorders, multimeters = {}, {}
        interval = max(0.1, self.cfg.dt)
        for layer, nodes in layer_nodes.items():
            if nodes is None or len(nodes) == 0:
                continue
            sr = nest.Create("spike_recorder")
            nest.Connect(nodes, sr)
            recorders[layer] = sr
            if self.sim.record_traces:
                # Record V_m from a sample of cells (up to 20) at >=0.5 ms so the
                # per-layer mean V_m forms a usable LFP proxy without huge logs.
                mm = nest.Create("multimeter", params={
                    "record_from": ["V_m"], "interval": max(0.5, self.cfg.dt)})
                nest.Connect(mm, nodes[: min(20, len(nodes))])
                multimeters[layer] = mm
        return recorders, multimeters

    @staticmethod
    def _events(device, nest):
        try:
            return device.get("events")
        except (AttributeError, TypeError):
            return nest.GetStatus(device, "events")[0]

    def _collect(self, recorders, multimeters):
        nest = self._nest
        spikes, traces = {}, {}
        for layer, sr in recorders.items():
            ev = self._events(sr, nest)
            order = np.argsort(np.asarray(ev["times"]))
            spikes[layer] = {
                "times": np.asarray(ev["times"])[order],
                "senders": np.asarray(ev["senders"])[order],
            }
        for layer, mm in multimeters.items():
            ev = self._events(mm, nest)
            times = np.asarray(ev["times"], float)
            volts = np.asarray(ev["V_m"], float)
            if not len(times):
                continue
            # mean V_m across the recorded cells at each sampled time -> LFP proxy
            ut = np.unique(times)
            idx = np.searchsorted(ut, times)
            sums = np.bincount(idx, weights=volts, minlength=len(ut))
            cnts = np.bincount(idx, minlength=len(ut))
            traces[layer] = {"time": ut, "voltage": sums / np.maximum(1, cnts)}
        return spikes, traces

    # -- run ---------------------------------------------------------------

    def run(self):
        """Build, drive, simulate. Returns (spikes, traces, meta)."""
        if not self._nest_available:
            raise RuntimeError(
                "NEST is not installed; this model requires NEST. "
                "Install NEST 3.x or run inside a NEST-enabled environment."
            )
        self._init_kernel()
        pops = self.build_network()
        self.connect_network(pops)
        self.attach_auditory_input(pops)
        self.attach_sleep_drives(pops)
        recorders, multimeters = self._setup_recorders(pops)

        self._nest.Simulate(self.cfg.tstop)

        spikes, traces = self._collect(recorders, multimeters)
        meta = {
            "tstop": self.cfg.tstop,
            "dt": self.cfg.dt,
            "neuron_model": self._neuron_model,
            "n_per_layer": {l: (len(n) if n is not None else 0)
                            for l, n in self._layer_nodes(pops).items()},
            "slow_freq": self.sleep.slow_freq,
            "spindle_freq": self.sleep.spindle_freq,
        }
        return spikes, traces, meta
