# Auditory thalamo-cortical sleep model (NEST)

A closed-loop thalamo-cortical column for the **auditory** system that, under
sleep drives, produces the **slow oscillation (~1 Hz)** with **sleep spindles
(~13 Hz)** nested on its UP states. Built directly in NEST, **no optimisation**.

The column architecture is adapted from the Cortical Column model
([max-talanov/cc, `optimization_nest`](https://github.com/max-talanov/cc/tree/main/optimization_nest),
`core/simulator.py`): five stages of `iaf_cond_exp` neurons — thalamus
(TCR excitatory / nRT inhibitory), L4, L2/3 (RS + FRB), L5 (TuftRS + TuftIB),
L6, with a shared L5/6 interneuron pool (Basket / LTS / Axoaxonic) — wired with
`pairwise_bernoulli` connections. The optimisation harness (`run_opt.py`,
`methods/`) is dropped; this builds one column and runs it.

## Neuron model

Every neuron (thalamic and cortical, excitatory and inhibitory) is NEST's
**`iaf_cond_exp`** — a leaky integrate-and-fire point neuron with
conductance-based, exponentially-decaying synapses; excitation vs. inhibition is
routed by the sign of the synaptic weight. The full governing equations, the
`aeif_cond_exp` fallback, and the parameter table (with Mushtaq 2024 sources)
are in **[`docs/model_equations.md`](docs/model_equations.md)**.

## What this adds over the cc column

1. **A closed thalamo-cortical loop.** cc is feed-forward (thalamus→L4, L6→L5).
   We add the two missing return projections that close the loop:
   - **corticothalamic feedback** `L6_E → TCR` and `L6_E → nRT`;
   - **reciprocal reticular inhibition** `nRT → TCR` and `nRT → nRT`.
   The TCR↔nRT reciprocal loop is the canonical spindle generator; its delays
   are tuned so the round trip resonates near 13 Hz.

2. **Sleep drives** (`SleepParams`):
   - **Slow oscillation (1 Hz):** an `ac_generator` injects a ~1 Hz current into
     the cortical pyramidal pools so the column cycles between depolarised **UP**
     (firing) and silent **DOWN** states.
   - **Spindles (13 Hz):** a ~13 Hz `ac_generator` into the thalamus. The
     thalamic DC is biased so relay/reticular cells are spindle-permissive
     (near threshold) only on the UP phase; the 13 Hz peaks then evoke
     UP-gated, waxing/waning spindles, reinforced by the nRT↔TCR loop.

The thalamus is read as the **MGB (medial geniculate body)**, the auditory
relay; its background drive (`poisson_generator`) represents auditory-nerve /
inferior-colliculus input (kept light, as the relay is rhythmic during sleep).

### Honest scope

`iaf_cond_exp` has no intrinsic Ih / T-type Ca²⁺ currents, so spindles here are
an **imposed-drive + loop-resonance hybrid**, not purely intrinsic thalamic
rhythmogenesis. Swapping the neuron model to NEST's `ht_neuron` (Hill–Tononi,
with Ih/IT/INaP/IKNa) would make both rhythms fully emergent — a documented
upgrade for future bio-plausible runs.

## Parameters from Mushtaq et al. 2024

The synaptic biophysics and network sizing are aligned with the conductance-based
thalamocortical sleep model of **Mushtaq, Marshall, ul Haq & Martinez (2024),
*"Possible mechanisms to improve sleep spindles via closed loop stimulation
during slow wave sleep: A computational study"*, PLOS ONE 19(6):e0306218**
([doi](https://doi.org/10.1371/journal.pone.0306218)). That model is full
Hodgkin–Huxley with PY/IN cortical and TC/RE thalamic cells, so its **absolute,
area-based conductances (mS/cm²·area) are not transferable** to our point
`iaf_cond_exp` cells. What *is* transferable was adopted:

- **Synaptic reversal potentials** (Table 3 / "Synaptic currents"): AMPA/NMDA
  `E_syn = 0 mV`, **GABA_A `E_syn = −70 mV`** (was −75). GABA_B / `E_K = −95 mV`
  is noted but lumped into the single GABA conductance (no separate GABA_B here).
- **Synaptic kinetics** from the first-order rate constants: **AMPA `τ ≈ 5.3 ms`**
  (β = 0.19 ms⁻¹; was 2.0) and **GABA_A `τ ≈ 6.0 ms`** (β = 0.166 ms⁻¹; already
  matched). NMDA (`τ ≈ 150 ms`) is not represented by the single-exponential
  conductance and is omitted. Because the slower AMPA τ alone drives the
  recurrent cortex into runaway, the intracortical excitatory weight is scaled
  by ≈ 2.0/5.3 (**charge-preserving calibration**, holding weight × τ roughly
  constant) so the network keeps the article kinetics in a plausible firing
  regime.
- **Relative thalamic-loop conductances** (Table 3, intra-thalamic block):
  TC→RE (AMPA) strongest, RE→TC (GABA_A + GABA_B) and RE→RE (GABA_A) — the
  resonator weights honour these ratios.
- **Corticothalamic feedback ratio** (Table 3, cortico-thalamic block):
  **PY→TC = 2 × PY→RE**, so L6→TCR feedback is twice L6→nRT.
- **TC→IN thalamocortical inhibition** (Table 3) — a feedforward TC→interneuron
  projection (more focal than TC→PY) that cc lacked, added as `thalamus_E → L4_I`.
- **Network sizing** (Fig 1/3, "Network geometry"): `config/network_auditory_mushtaq.yaml`
  uses **PY = 200 / IN = 40 / TC = 40 / RE = 40**, with the 200 excitatory cells
  distributed across the auditory laminae (E:I = 5:1 preserved) and the thalamus
  matching exactly. At this sizing the model reproduces the article's control
  rhythms cleanly: **slow wave ≈ 1.00 Hz, spindle ≈ 13.0 Hz** (within Mushtaq's
  10–16 Hz control band).

## Files

| file | purpose |
|------|---------|
| `tc_network.py` | `AuditoryThalamoCorticalSleep` — builds/wires/drives/runs the column |
| `tc_run.py` | CLI driver: run, build LFP-proxy signals, verify 1 Hz & 13 Hz, plot |
| `tc_architecture.py` | draws the thalamo-cortical loop architecture schematic (`out/tc_architecture.png`, no NEST needed) |
| `config/network_auditory_local.yaml` | small (~112 neurons), short — local sane test |
| `config/network_auditory_mushtaq.yaml` | Mushtaq et al. 2024 sizing (PY 200 / IN 40 / TC 40 / RE 40), 15 s |
| `config/network_auditory_mn5.yaml` | realistic (~1150 neurons), long — bio-plausible MN5 |
| `slurm/tc_sleep_mn5.sbatch` | MareNostrum 5 submission (single long, multi-threaded run) |

## Run

Local sane test (seconds on a laptop with NEST 3.x installed):

```bash
pip install -r ../requirements.txt   # numpy matplotlib scipy pyyaml (+ NEST)
python3 tc_sleep/tc_run.py --config tc_sleep/config/network_auditory_local.yaml --outdir out
```

It prints per-layer firing rates and the **detected slow-wave (~1 Hz)** and
**spindle (~13 Hz)** peaks, and writes `out/tc_sleep_local.png` with five
indexed panels:
- **(a)** per-layer spike raster with UP/DOWN banding;
- **(b)** an **EEG/LFP-like composite trace built from the recorded membrane
  potentials** — cortical mean V_m (low-pass <2 Hz) for the slow wave plus
  thalamic mean V_m (9–16 Hz band-pass) for the ~13 Hz spindle wavelets,
  superimposed the way a sleep EEG channel looks — with detected **SP**
  (spindle) and **SW** (slow-wave) epochs shaded;
- **(c)** a zoom over ~2 slow cycles resolving the individual 13 Hz spindle
  wavelets riding on the slow wave (slow V_m component overlaid in black);
- **(d)** a thalamic **V_m** spectrogram (13 Hz bursts gated to UP states);
- **(e)** cortical and thalamic **V_m** PSDs (slow and spindle peaks).

The per-layer mean V_m used for the figure is an LFP proxy taken from
multimeters on a sample of cells; the printed/asserted slow-wave and spindle
peaks are computed independently from the population spike rates, so the
self-validation still reflects that neurons actually fire at those rhythms.

It also writes `out/tc_sleep_<tag>_decomp.png`, an **integrated single-graph
decomposition** in the style of intracranial sleep-LFP figures (e.g. a Reuniens
LFP trace, or the SP/SW EEG panels): one composite auditory-cortex LFP on top,
then the **Spindle (7–15 Hz)** and **Slow wave (0.5–2 Hz)** bands band-passed
*out of that same trace* and stacked beneath it, each with a µV/300 ms scale
bar and the detected **SP**/**SW** epochs shaded. Because both bands are
extracted from one signal, the figure shows directly that the two rhythms
coexist in a single channel, with the 7–15 Hz spindles waxing and waning on the
slow-wave UP states (`make_decomposition_plot` in `tc_run.py`).

The process exits non-zero if either rhythm is missing, so the run
self-validates.

MareNostrum 5 bio-plausible run:

```bash
sbatch tc_sleep/slurm/tc_sleep_mn5.sbatch     # edit account/qos/partition + module loads first
```

## Scaling (local ↔ MN5)

With fixed `pairwise_bernoulli` probability, per-neuron in-degree grows with
population size, so a 10× larger column would saturate the cortex. Synaptic
weights are therefore scaled by `REF_E / N_L23_E` (balanced-network scaling),
keeping per-neuron input and the E/I balance scale-invariant. The factor is
`1.0` for the local default sizes, so the local results are unaffected, while
the MN5 column reproduces the same firing regime at ~10× the neurons.
