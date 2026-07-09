# Neuron model and governing equations

The simulation supports **two selectable neuron models** (`--neuron-model`, or
`simulation.neuron_model` in a config), used uniformly across all populations —
thalamic (TCR/relay, nRT/reticular) and cortical (L4, L2/3, L5, L6):

1. **`iaf_cond_exp`** (default) — a **leaky integrate-and-fire (LIF)** point
   neuron with conductance-based, exponential synapses. Rhythms are an
   **imposed-drive + loop-resonance hybrid**.
2. **`ht_neuron`** (Hill & Tononi 2005) — a **Hodgkin–Huxley** neuron with
   intrinsic currents, giving **emergent** rhythms. See
   [§ ht_neuron](#ht_neuron-hodgkinhuxley-emergent-rhythms) below.

The default preference list in [`tc_network.py`](../tc_network.py):

```python
_PREFERRED_NEURON_MODELS = ["iaf_cond_exp", "aeif_cond_exp"]
```

`iaf_cond_exp` routes excitation vs. inhibition by the **sign of the synaptic
weight** (positive → `g_ex`/`E_ex`, negative → `g_in`/`E_in`); `ht_neuron`
routes by **receptor type** (AMPA/NMDA/GABA_A/GABA_B).

> **Scope note (`iaf_cond_exp`).** The default point LIF neuron has **no
> intrinsic Iₕ, T-type Ca²⁺, or persistent-Na currents**, so under it the slow
> oscillation and spindles are imposed/loop-resonant rather than emergent. The
> **`ht_neuron`** model below adds exactly those currents and makes both rhythms
> emergent — the key modelling upgrade toward Mushtaq et al. 2024's full
> Hodgkin–Huxley cells.

## `iaf_cond_exp` equations

### Subthreshold membrane dynamics

$$
C_m \frac{dV}{dt} = -g_L\,\bigl(V - E_L\bigr)
        - g_{\mathrm{ex}}(t)\,\bigl(V - E_{\mathrm{ex}}\bigr)
        - g_{\mathrm{in}}(t)\,\bigl(V - E_{\mathrm{in}}\bigr)
        + I_e(t)
$$

### Conductance-based exponential synapses

$$
\frac{dg_{\mathrm{ex}}}{dt} = -\frac{g_{\mathrm{ex}}}{\tau_{\mathrm{syn,ex}}},
\qquad
\frac{dg_{\mathrm{in}}}{dt} = -\frac{g_{\mathrm{in}}}{\tau_{\mathrm{syn,in}}}
$$

### Synaptic jump (arrival of presynaptic spike $k$ with weight $w_k$)

$$
g_{\mathrm{ex}} \leftarrow g_{\mathrm{ex}} + w_k \quad (w_k > 0),
\qquad
g_{\mathrm{in}} \leftarrow g_{\mathrm{in}} + |w_k| \quad (w_k < 0)
$$

### Threshold and reset

$$
\text{if } V(t^-) \ge V_{\mathrm{th}}:\quad
   V \leftarrow V_{\mathrm{reset}}, \quad
   V(t)=V_{\mathrm{reset}}\ \text{for } t\in[t^*,\,t^*+t_{\mathrm{ref}}]
$$

### Injected drive current (slow oscillation + spindle AC generators)

$$
I_e(t) = I_{\mathrm{slow}} + A_{\mathrm{slow}}\sin\!\bigl(2\pi f_{\mathrm{slow}} t\bigr)
       + A_{\mathrm{spin}}\sin\!\bigl(2\pi f_{\mathrm{spin}} t + \varphi\bigr)
$$

The auditory background is delivered as an additional excitatory conductance
driven by a `poisson_generator` onto the relay cells (same $g_{\mathrm{ex}}$
pathway as above), representing auditory-nerve / inferior-colliculus input.

## Fallback: `aeif_cond_exp` (adaptive exponential I&F)

Used only if `iaf_cond_exp` is unavailable. It adds a spike-generating
exponential term and an adaptation current $w$:

$$
C_m \frac{dV}{dt} = -g_L(V-E_L)
    + g_L\,\Delta_T\,\exp\!\Bigl(\tfrac{V-V_T}{\Delta_T}\Bigr)
    - g_{\mathrm{ex}}(t)(V-E_{\mathrm{ex}})
    - g_{\mathrm{in}}(t)(V-E_{\mathrm{in}})
    - w + I_e(t)
$$

$$
\tau_w \frac{dw}{dt} = a\,(V-E_L) - w,
\qquad
w \leftarrow w + b \ \text{at each spike}
$$

## `ht_neuron` (Hodgkin–Huxley, emergent rhythms)

Selected with `--neuron-model ht_neuron` (or `config/network_auditory_hh.yaml`).
This is NEST's implementation of the **Hill & Tononi (2005)** thalamocortical
neuron. The membrane potential integrates conductance-based synaptic currents,
a spike-generating term, and a set of **intrinsic currents**:

$$
\tau_m \frac{dV}{dt} = -\bigl(V - E_L\bigr)
   - \frac{1}{g_{\mathrm{NaL}}+g_{\mathrm{KL}}}\Bigl(
       I_{\mathrm{Na}} + I_{\mathrm{K}}
     + I_{\mathrm{NaP}} + I_{\mathrm{KNa}}
     + I_{T} + I_{h}
     + I_{\mathrm{syn}} \Bigr)
$$

with the intrinsic currents (each $I_x = g_x\,m^{M}h^{N}(V-E_x)$):

| current | symbol | role | present in |
|---------|--------|------|------------|
| persistent Na⁺ | $I_{\mathrm{NaP}}$ | depolarising drive for UP states | cortex (PY/IN) |
| Na-dependent K⁺ | $I_{\mathrm{KNa}}$ | spike-frequency adaptation → ends UP states | cortex (PY/IN) |
| low-threshold Ca²⁺ | $I_{T}$ | post-inhibitory **rebound burst** (the spindle) | TC, RE |
| h-current | $I_{h}$ | pacemaker; waxing/waning + termination | TC (relay) |

The **spindle** is emergent: reticular (RE) GABA_A/GABA_B inhibition
hyperpolarises the relay (TC) cell, de-inactivating $I_T$; on release, the $I_T$
rebound fires TC, which re-excites RE (AMPA), and the loop rings in the
spindle band. The **slow oscillation** is emergent from cortical
$I_{\mathrm{NaP}}$ (depolarising) balanced by $I_{\mathrm{KNa}}$ adaptation.

Synaptic currents are receptor-typed beta-functions with reversal potentials
$E_{\mathrm{AMPA}}=E_{\mathrm{NMDA}}=0$, $E_{\mathrm{GABA_A}}=-70$,
$E_{\mathrm{GABA_B}}=-90$ mV:

$$
I_{\mathrm{syn}} = \sum_{r\in\{\mathrm{AMPA,NMDA,GABA_A,GABA_B}\}}
    g_{\mathrm{peak},r}\,w\,s_r(t)\,\bigl(V - E_r\bigr)
$$

The sleep/wake regime is set by the K-leak conductance $g_{\mathrm{KL}}$
(neuromodulation): thalamic cells use a burst-permissive value. No 13 Hz drive
is injected in this mode; the light 1 Hz `ac_generator` only gates the emergent
spindles to UP states. The intrinsic conductances and loop calibration live in
the `HHParams` dataclass in [`tc_network.py`](../tc_network.py).

## Synaptic neurotransmission

Both neuron models use **conductance-based** synapses: a presynaptic spike opens
postsynaptic receptor channels, and the resulting current depends on the driving
force $(V - E_r)$ toward that receptor's reversal potential $E_r$. The total
synaptic current entering a cell is the sum over all its receptor types $r$:

$$
I_{\mathrm{syn}}(t) = \sum_{r} g_r(t)\,f_r(V)\,\bigl(V - E_r\bigr)
$$

where $g_r(t)$ is the (time-varying) open conductance of receptor $r$, $E_r$ its
reversal potential (AMPA/NMDA $= 0$, GABA$_A = -70$, GABA$_B/E_K = -90{\dots}{-}95$
mV), and $f_r(V)$ a voltage-gating factor ($f_r \equiv 1$ except for NMDA). The
two models differ only in **how $g_r(t)$ is built from the incoming spikes**.

### 1. `iaf_cond_exp` — single-exponential conductances

Two lumped receptor conductances, excitatory ($g_{\mathrm{ex}}$, "AMPA-like") and
inhibitory ($g_{\mathrm{in}}$, "GABA$_A$-like"). Each presynaptic spike $k$ of
weight $w_k$ arriving at time $t_k$ **instantaneously steps** the conductance,
which then **decays exponentially**:

$$
g_{\mathrm{ex}}(t) = \sum_{k}\, w_k \;\exp\!\Bigl(-\tfrac{t - t_k}{\tau_{\mathrm{syn,ex}}}\Bigr)\,\Theta(t-t_k),
\qquad
\frac{dg_{\mathrm{ex}}}{dt} = -\frac{g_{\mathrm{ex}}}{\tau_{\mathrm{syn,ex}}}
$$

(identically for $g_{\mathrm{in}}$ with $\tau_{\mathrm{syn,in}}$), $\Theta$ the
Heaviside step. Excitation vs. inhibition is chosen by the **sign of $w_k$**
(positive $\to g_{\mathrm{ex}}$, negative $\to g_{\mathrm{in}}$). The postsynaptic
current is then $I_{\mathrm{syn}} = g_{\mathrm{ex}}(V-E_{\mathrm{ex}}) + g_{\mathrm{in}}(V-E_{\mathrm{in}})$.
Time constants $\tau_{\mathrm{syn,ex}}=5.3$ ms and $\tau_{\mathrm{syn,in}}=6.0$ ms
(Mushtaq 2024; see below).

### 2. `ht_neuron` — beta-function conductances (rise + decay)

Four receptor types are modelled separately — **AMPA, NMDA, GABA$_A$, GABA$_B$**
— each with its own **rise and decay** kinetics (a beta function / difference of
two exponentials), routed by `receptor_type` rather than weight sign. A single
presynaptic spike of weight $w$ produces

$$
g_r(t) = w\,\bar{g}_r\,B_r(t),
\qquad
B_r(t) = \frac{\exp\!\bigl(-t/\tau_{\mathrm{decay},r}\bigr) - \exp\!\bigl(-t/\tau_{\mathrm{rise},r}\bigr)}{\;\exp\!\bigl(-t_{\mathrm{peak},r}/\tau_{\mathrm{decay},r}\bigr) - \exp\!\bigl(-t_{\mathrm{peak},r}/\tau_{\mathrm{rise},r}\bigr)\;}
$$

normalised so that $\max_t B_r(t) = 1$ (hence $\bar{g}_r = g_{\mathrm{peak},r}$ is
the **peak** conductance per unit weight), with peak time

$$
t_{\mathrm{peak},r} = \frac{\tau_{\mathrm{rise},r}\,\tau_{\mathrm{decay},r}}{\tau_{\mathrm{decay},r}-\tau_{\mathrm{rise},r}}\,\ln\!\frac{\tau_{\mathrm{decay},r}}{\tau_{\mathrm{rise},r}}.
$$

Conductances **sum linearly** over all incoming spikes (each shifted to its own
arrival time). Default kinetics: AMPA $\tau_{\mathrm{rise}}/\tau_{\mathrm{decay}}
= 0.5/2.4$ ms, GABA$_A = 1.0/7.0$ ms, NMDA $= 4.0/40$ ms, and the **slow**
GABA$_B = 60/200$ ms — the long GABA$_B$ tail is what paces the spindle
waxing/waning on RE→TC.

**NMDA voltage gating (Mg²⁺ block).** NMDA carries an extra voltage-dependent
factor $f_{\mathrm{NMDA}}(V) = m(V)$ that relieves the Mg²⁺ block as the cell
depolarises:

$$
m_\infty(V) = \frac{1}{1 + \exp\!\bigl(-S_{\mathrm{act}}\,(V - V_{\mathrm{act}})\bigr)},
\qquad S_{\mathrm{act}} = 0.081\ \mathrm{mV^{-1}},\ V_{\mathrm{act}} = -25.6\ \mathrm{mV}.
$$

With `instant_unblock_NMDA = True` (as used here) $m(V)=m_\infty(V)$ follows the
voltage instantaneously; otherwise $m$ relaxes toward $m_\infty$ with fast/slow
time constants $\tau_{\mathrm{Mg,fast}}=0.68$, $\tau_{\mathrm{Mg,slow}}=22.7$ ms.
For all other receptors $f_r(V)\equiv 1$.

Putting it together, the `ht_neuron` synaptic current is

$$
I_{\mathrm{syn}} = g_{\mathrm{AMPA}}(V-E_{\mathrm{AMPA}})
   + m(V)\,g_{\mathrm{NMDA}}(V-E_{\mathrm{NMDA}})
   + g_{\mathrm{GABA_A}}(V-E_{\mathrm{GABA_A}})
   + g_{\mathrm{GABA_B}}(V-E_{\mathrm{GABA_B}}).
$$

> **Which receptors this model actually wires.** Excitatory projections use
> **AMPA**; inhibition uses **GABA$_A$**, with an added slow **GABA$_B$**
> component on the reticular→relay (RE→TC) synapse (Mushtaq 2024, Table 3).
> NMDA is available in `ht_neuron` but is not currently wired. The per-synapse
> weight $w$ scales $g_{\mathrm{peak},r}$; the intra-thalamic loop weights are
> **fan-in normalised** so the total conductance onto each cell is
> size-invariant (see `_ht_thal_w` in [`tc_network.py`](../tc_network.py)).

### 3. Biophysical reference (Mushtaq et al. 2024)

The reversal potentials and decay time constants above come from Mushtaq 2024,
whose conductance-based model computes the open fraction $[O]$ with a
**first-order kinetic scheme** driven by a brief transmitter pulse $[T]$:

$$
I_{\mathrm{syn}} = g_{\mathrm{syn}}\,[O]\,f(V)\,(V - E_{\mathrm{syn}}),
\qquad
\frac{d[O]}{dt} = \alpha\,(1-[O])\,[T] - \beta\,[O],
$$

$$
[T] = A\,\Theta(t_0 + t_{\max} - t)\,\Theta(t - t_0),
\qquad t_{\max}=0.3\ \mathrm{ms},\ A=0.5,
$$

with rate constants (ms$^{-1}$): AMPA $\alpha{=}1.1,\ \beta{=}0.19$; GABA$_A$
$\alpha{=}10.5,\ \beta{=}0.166$; NMDA $\alpha{=}1,\ \beta{=}0.0067$. The decay
time constant is $\tau_{\mathrm{decay}} = 1/\beta$ — this is where our
$\tau_{\mathrm{syn,ex}}\approx 5.3$ ms and $\tau_{\mathrm{syn,in}}\approx 6.0$ ms
originate. Intracortical AMPA/GABA$_A$ additionally carry short-term depression
$D_{n+1} = 1 - (1 - D_n(1-U))\,e^{-\Delta t/\tau}$ ($\tau=700$ ms, $U=0.07/0.073$),
and **GABA$_B$** uses a higher-order G-protein scheme,

$$
I_{\mathrm{GABA_B}} = g_{\mathrm{GABA_B}}\,\frac{[G]^4}{[G]^4 + K_d}\,(V - E_K),
\qquad
\frac{d[R]}{dt} = K_1(1-[R])[T] - K_2[R],
\qquad
\frac{d[G]}{dt} = K_3[R] - K_4[G],
$$

($E_K=-95$ mV; $K_1{=}0.052,\ K_2{=}0.0013,\ K_3{=}0.098,\ K_4{=}100$). Because
our cells are point/`ht_neuron` rather than Mushtaq's multi-compartment HH, these
schemes are represented by the exponential (`iaf_cond_exp`) and beta-function
(`ht_neuron`) conductances above rather than integrated verbatim, but the
reversals, kinetics, and relative conductances are carried over (see the README's
"Parameters from Mushtaq et al. 2024").

## Parameter values

Set in [`_neuron_params()`](../tc_network.py) and `SynapseParams`
(post-Mushtaq-2024 update). The membrane time constant is
$\tau_m = C_m/g_L = 100\,\mathrm{pF} / 10\,\mathrm{nS} = 10\ \mathrm{ms}$.

| Symbol | Meaning | Value | Source |
|--------|---------|-------|--------|
| $C_m$ | membrane capacitance | 100 pF | cc |
| $g_L$ | leak conductance | 10 nS | cc |
| $E_L$ | leak reversal | −70 mV | cc |
| $V_{\mathrm{th}}$ | spike threshold | −55 mV | cc |
| $V_{\mathrm{reset}}$ | reset potential | −60 mV | cc |
| $t_{\mathrm{ref}}$ | refractory period | 2 ms | cc |
| $E_{\mathrm{ex}}$ | excitatory (AMPA) reversal | 0 mV | Mushtaq 2024, Table 3 |
| $E_{\mathrm{in}}$ | inhibitory (GABA_A) reversal | −70 mV | Mushtaq 2024, Table 3 |
| $\tau_{\mathrm{syn,ex}}$ | AMPA decay | 5.3 ms | Mushtaq 2024 (β = 0.19 ms⁻¹) |
| $\tau_{\mathrm{syn,in}}$ | GABA_A decay | 6.0 ms | Mushtaq 2024 (β = 0.166 ms⁻¹) |
| $f_{\mathrm{slow}}$ | slow-oscillation drive | 1 Hz | — |
| $f_{\mathrm{spin}}$ | spindle drive | 13 Hz | — |

## References

- **`iaf_cond_exp`** — NEST integrate-and-fire, conductance-based, exponential
  synaptic currents. See the [NEST model documentation](https://nest-simulator.readthedocs.io/en/stable/models/iaf_cond_exp.html).
- **Mushtaq, Marshall, ul Haq & Martinez (2024)**, *Possible mechanisms to
  improve sleep spindles via closed loop stimulation during slow wave sleep: A
  computational study*, PLOS ONE 19(6):e0306218.
  [doi:10.1371/journal.pone.0306218](https://doi.org/10.1371/journal.pone.0306218)
  — source of the synaptic reversal potentials, kinetics, and network sizing
  (see the "Parameters from Mushtaq et al. 2024" section of the README).
