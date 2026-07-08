# Neuron model and governing equations

All neurons in this simulation — thalamic (TCR relay, nRT reticular) and cortical
(L4, L2/3, L5, L6; excitatory and inhibitory) — are the **same** point neuron
model: NEST's **`iaf_cond_exp`**, a **leaky integrate-and-fire (LIF) neuron with
conductance-based, exponentially-decaying synapses**. Subtypes differ only in
their connectivity and drive, not in their intrinsic dynamics.

The model is selected in [`tc_network.py`](../tc_network.py):

```python
_PREFERRED_NEURON_MODELS = ["iaf_cond_exp", "aeif_cond_exp"]
```

`iaf_cond_exp` is used when available; `aeif_cond_exp` (adaptive exponential
integrate-and-fire) is only a fallback. Excitation vs. inhibition is routed by
the **sign of the synaptic weight** (positive → `g_ex`/`E_ex`, negative →
`g_in`/`E_in`), not by a receptor type.

> **Scope caveat.** This is a *point* LIF neuron: it has **no intrinsic Iₕ,
> T-type Ca²⁺, or persistent-Na currents**. This is the key modelling difference
> from Mushtaq et al. 2024, whose cells are full Hodgkin–Huxley. Here the slow
> oscillation and spindles are an **imposed-drive + loop-resonance hybrid**
> rather than emergent from intrinsic currents. Swapping the neuron model to
> NEST's `ht_neuron` (Hill–Tononi) would make both rhythms emergent.

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
