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

## Files

| file | purpose |
|------|---------|
| `tc_network.py` | `AuditoryThalamoCorticalSleep` — builds/wires/drives/runs the column |
| `tc_run.py` | CLI driver: run, build LFP-proxy signals, verify 1 Hz & 13 Hz, plot |
| `tc_architecture.py` | draws the thalamo-cortical loop architecture schematic (`out/tc_architecture.png`, no NEST needed) |
| `config/network_auditory_local.yaml` | small (~112 neurons), short — local sane test |
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
