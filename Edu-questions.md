# Educational questions about the project

## 1) How do we know when a spindle starts? (grey rectangles in the figure)

Figure: `out/compare_regularity_vs_literature_score.png` (commit `7b3e1c8`).

The grey rectangles come from a simple threshold on the spindle-band amplitude.
The figure is drawn by `make_comparison_figure` in `neuron/ctx_analyze.py`
(committed in `7170e23`), which detects events with `_detect_spindles`. The
same definition is used for scoring in `_literature_score` in
`neuron/ctx_thalamus_mpi.py`; the two are kept in sync by hand.

### How a spindle start is found

1. **Build an LFP-like signal.** For each excitatory layer (L2/3, L4, L5, L6),
   spike times are binned at 1 ms. Each binned signal is convolved with a
   synaptic-shaped kernel (2 ms rise, 100 ms decay) and z-scored. L5 and L6 are
   sign-flipped, and the four layers are averaged into one signal.
2. **Band-pass the signal to 10–15 Hz.** This is the lower trace of each pair
   in the figure.
3. **Take the Hilbert envelope**, `env = |hilbert(spin)|`, which is the
   instantaneous amplitude of that oscillation.
4. **Set the threshold** at `mean(env) + 1.5·SD(env)`, computed over the whole
   run except the first second. That second holds the network's start-up
   transient and filter edge effects, and it is excluded from detection too.
   The exclusion was added after the comparison figure was made, so the
   committed PNG and earlier grid-sweep numbers used the whole run.
5. **Mark events:** a spindle starts where the envelope rises above the
   threshold and ends where it falls back below. Each grey rectangle is one
   span where the envelope stays above the threshold.

### What this means for the figure

- **The "start" is when the envelope crosses the threshold, not when the
  spindle visibly begins.** In Option 1 the oscillation can be seen building up
  before the grey box. The box only covers the high-amplitude peak, roughly
  100–150 ms.
- **There is no minimum duration and no merging of nearby events.** That's why
  a very thin sliver appears near 9.6 s in Option 2. Two short crossings close
  together count as two events.
- **The threshold adapts to each run.** Option 2 has 10–15 Hz activity almost
  everywhere, which raises the envelope's mean and SD. As a result, only the
  few biggest bursts cross the threshold. So the low event count in Option 2
  partly reflects this choice of threshold, not only how the network behaves.
- **This is simpler than standard spindle detection in the literature.**
  Detectors there usually require a duration of about 0.5–2 s and use two
  thresholds: a high one to detect the event and a lower one to set its
  boundaries. The events here are much shorter than that. So the shaded boxes
  are really "spindle-band amplitude peaks", not spindles in the strict EEG
  sense. That matters if `n_spindle_events` or `spindle_cv` are compared
  directly to published numbers.

### Possible follow-ups

- Add two-threshold detection with a minimum duration, so start and end times
  match the literature definition.

## 2) What is the network architecture? (number of neurons and synapses/projections)

Source of truth: `ParallelCorticoThalamicNet` in `neuron/ctx_thalamus_mpi.py`
(the MPI/MN5 model). All counts below are computed from its wiring code
(`DEFAULT_SIZES`, `_project`, `_wire_re_re_local`, `_wire_l5_recurrent`,
`_wire_gap_range`) with the production defaults `--conv 100`, `ib_frac=0.5`,
`gap_deg=6`, `gap_short=2`.

Two sizes are in use:

- `--scale 1.0` gives **3050 cells**, the sizes declared in
  `config/network_auditory_mn5_5k.yaml`.
- `--scale 1.65` gives **5031 cells**. This is the MN5 full-scale run (job
  45453403) and the size shown in the diagram.

### Diagram (scale 1.65, 5031 cells)

Solid arrows are excitatory (AMPA-like, E = 0 mV). Dotted arrows are inhibitory
(GABA_A). Double-headed thick links are gap junctions. Edge labels give the
number of connections (NetCons) in the projection.

```mermaid
flowchart TB
    subgraph THAL["Thalamus (437 cells)"]
        direction LR
        TC["TC relay<br/>346"]
        RE["RE reticular<br/>91"]
    end

    subgraph CTX["Cortical column (4594 cells)"]
        direction TB
        subgraph L23["L2/3"]
            direction LR
            L23E["L2/3 E<br/>1056"]
            L23I["L2/3 I<br/>264"]
        end
        subgraph L4["L4"]
            direction LR
            L4E["L4 E<br/>1056"]
            L4I["L4 I<br/>264"]
        end
        subgraph L5["L5"]
            direction LR
            L5E["L5 E<br/>874<br/>(437 IB + 437 RS)"]
            L5I["L5 I<br/>214"]
        end
        subgraph L6["L6"]
            direction LR
            L6E["L6 E<br/>693"]
            L6I["L6 I<br/>173"]
        end
    end

    %% thalamic loop
    TC -->|"9.1k"| RE
    RE -.->|"31.5k (all-to-all, k=91)"| TC
    RE -.->|"358 (±2 neighbours)"| RE
    RE <==>|"gap 715"| RE

    %% ascending (thalamocortical)
    TC -->|"105.6k"| L4E

    %% feedforward intracortical
    L4E -->|"105.6k"| L23E
    L4E -->|"87.4k"| L5E
    L23E -->|"87.4k"| L5E
    L5E -->|"105.6k"| L23E
    L5E -->|"69.3k"| L6E

    %% L5 recurrence (Option 1/2, opt-in) + IB gap junctions
    L5E -->|"87.4k NMDA (opt-in)"| L5E
    L5E <==>|"gap 3481 (IB only)"| L5E

    %% local E<->I per layer
    L4E -->|"26.4k"| L4I
    L4I -.->|"105.6k"| L4E
    L23E -->|"26.4k"| L23I
    L23I -.->|"105.6k"| L23E
    L5E -->|"21.4k"| L5I
    L5I -.->|"87.4k"| L5E
    L6E -->|"17.3k"| L6I
    L6I -.->|"69.3k"| L6E

    %% descending (corticothalamic)
    L6E -->|"34.6k"| TC
    L6E -->|"9.1k"| RE

    classDef exc fill:#e6f4ea,stroke:#1e8e3e,color:#000
    classDef inh fill:#fce8e6,stroke:#c5221f,color:#000
    class TC,L23E,L4E,L5E,L6E exc
    class RE,L23I,L4I,L5I,L6I inh
```

### Neurons

| Population | Cell model | scale 1.0 | scale 1.65 |
|---|---|---:|---:|
| TC (thalamocortical relay) | `tc_neuron.TCCell` | 210 | 346 |
| RE (thalamic reticular) | `tc_neuron.RECell` | 55 | 91 |
| L4 E / I | `PYCell` / `FSCell` | 640 / 160 | 1056 / 264 |
| L2/3 E / I | `PYCell` / `FSCell` | 640 / 160 | 1056 / 264 |
| L5 E / I | half `PYCellIB` (intrinsically bursting), half `PYCell` (RS) / `FSCell` | 530 / 130 | 874 / 214 |
| L6 E / I | `PYCell` / `FSCell` | 420 / 105 | 693 / 173 |
| **Total** | | **3050** | **5031** |

### Projections (chemical synapses)

Wiring uses **fixed convergence**. Each target cell draws `k = min(conv, N_pre)`
random presynaptic cells, so a projection has `N_post × k` connections. Each
target has one `Exp2Syn` per projection, and every connection is a NetCon onto
it with weight `g / k`. So `g` is the total conductance a target cell gets from
that projection, in µS. The shared delay is 1 ms.

| Projection | Type (E_rev, τ1/τ2 ms) | g (µS) | k | # conn. scale 1.0 | # conn. scale 1.65 |
|---|---|---:|---:|---:|---:|
| TC → RE | AMPA (0, 0.5/2) | 0.011 | 100 | 5 500 | 9 100 |
| RE → TC | GABA_A (−85, 1/8) | 0.015 | all RE (55 / 91) | 11 550 | 31 486 |
| RE → RE | GABA_A (−75, 1/6), ±2 neighbours | 0.006 ± 0.002 per conn. | ≤4 | 214 | 358 |
| TC → L4E | AMPA (0, 0.5/2) | 0.02 | 100 | 64 000 | 105 600 |
| L4E → L2/3E | AMPA | 0.0015 | 100 | 64 000 | 105 600 |
| L4E → L5E | AMPA | 0.00075 | 100 | 53 000 | 87 400 |
| L2/3E → L5E | AMPA | 0.0015 | 100 | 53 000 | 87 400 |
| L5E → L2/3E | AMPA | 0.0015 | 100 | 64 000 | 105 600 |
| L5E → L6E | AMPA | 0.009 | 100 | 42 000 | 69 300 |
| L6E → TC | AMPA | 0.03 | 100 | 21 000 | 34 600 |
| L6E → RE | AMPA | 0.03 | 100 | 5 500 | 9 100 |
| E → I, each layer (L4, L2/3, L5, L6) | AMPA | 0.02 | 100 | 16 000 / 16 000 / 13 000 / 10 500 | 26 400 / 26 400 / 21 400 / 17 300 |
| I → E, L4, L2/3, L6 | GABA_A (−75, 1/8) | 0.08 | 100 | 64 000 / 64 000 / 42 000 | 105 600 / 105 600 / 69 300 |
| I → E, L5 | GABA_A (−75, 1/8) | `g_i_e_l5`: 0 by default, 0.02 in Options 1/2 | 100 | 53 000 | 87 400 |
| L5E → L5E recurrence (opt-in) | NMDA (Options 1/2) or Exp2Syn | `g_l5_rec`: 0 by default, 0.013 in Options 1/2 | 100 | 53 000 | 87 400 |
| **Total without L5 recurrence** | | | | **662 264** | **1 104 944** |
| **Total with L5 recurrence** | | | | **715 264** | **1 192 344** |

The L5 I → E connections are always created, even when `g_i_e_l5 = 0`. They
just carry zero weight.

Options 1 and 2 are the two L5-recurrence tunings compared in question 1.
Both use `g_l5_rec = 0.013` with the NMDA mechanism. Option 1 uses
`tau2 = 115, taur = 100`; Option 2 uses `tau2 = 200, taur = 120`.

### Gap junctions (electrical)

Small-world topology: a ring with ±3 neighbours plus 2 random shortcuts per
cell. The counts below are one-directional `GapMPI` halves.

| Group | g (µS) | scale 1.0 | scale 1.65 |
|---|---:|---:|---:|
| RE ↔ RE | 0.03 | 423 | 715 |
| L5 IB ↔ L5 IB (the first half of L5E) | 0.02 | 2 103 | 3 481 |

### What the architecture does *not* contain

- Within-layer E → E synapses, except the opt-in L5 recurrence.
- Direct L2/3 → thalamus or L6 → L4 connections. The only cortical output to
  the thalamus is L6E.
- Any input to L4E other than TC. That makes TC → L4E the only ascending arm
  of the loop.
- Any interneuron subtypes. All I cells are a single fast-spiking type
  (`FSCell`). The NEST reference model splits cells into Basket/LTS/Axoaxonic
  and RS/FRB/Tuft subtypes, and this model does not reproduce those splits.

## 3) How is the LFP calculated from the neuronal activity?

The simulation doesn't compute a real LFP, because it never records any
voltages or currents. The only output saved is spike times:
`pc.spike_record` in `neuron/ctx_thalamus_mpi.py` (line ~188), written to the
`.npz` with gids and population ranges. The "LFP" is a **proxy built
afterwards from those spikes**. The code is `_synaptic_lfp` in
`neuron/ctx_analyze.py` (line ~242), with an identical copy
(`_synaptic_lfp_local`) in `neuron/ctx_thalamus_mpi.py` (line ~603).

### Steps

1. **One trace per layer.** The excitatory cells of each layer (L2/3E, L4E,
   L5E, L6E) are taken separately. Inhibitory cells are left out.
2. **Count spikes per 1 ms bin.** This gives the whole population's spike count
   over time, not per-cell rates.
3. **Smooth with a synaptic-shaped kernel.** Each spike is replaced by a
   PSP-like waveform, `exp(-t/100) − exp(-t/2)`: 2 ms rise, 100 ms decay, cut
   off at 800 ms and scaled to a peak of 1. The smoothing only uses past
   spikes. The idea is that a real LFP comes mostly from slow synaptic
   currents, not from the spikes themselves.
4. **Normalise each layer** to mean 0 and SD 1.
5. **Flip the sign of L5 and L6.** This mimics the usual dipole in laminar
   recordings, negative near the surface and positive in deep layers. In the
   LFP figure the RE trace is flipped too.
6. **Average the four layers** into one cortical signal. The thalamic trace is
   the average of TC and RE.
7. **Band-pass filter** (3rd-order Butterworth, applied forwards and backwards
   so there's no phase shift):
   - slow oscillation: 0.5–2 Hz
   - spindle band: 8–15 Hz in `make_lfp_figure`, but 10–15 Hz in
     `_literature_score` and in the comparison figure from question 1
   - Hilbert envelope of the spindle band

The "raw (SO+spindle)" panel in the comparison figure is **not** the unfiltered
signal: it is the 0.5–2 Hz band plus the 10–15 Hz band added together
(`make_literature_reconstruction_figure` / `make_comparison_figure` in
`neuron/ctx_analyze.py`). Because everything outside those two bands is
removed, any run will look like slow oscillations with spindles on top, so its
visual similarity to published raw traces is weak evidence on its own. The
code now labels this panel "SO + spindle bands (sum)"; the committed PNGs
made before that change still show the old "raw" label.

The mean-field figure uses a simpler proxy still: the population firing rate,
smoothed with a 5 ms moving average (`make_meanfield_figure` in
`neuron/ctx_analyze.py`).

### Example: LFP proxy from two full-scale MN5 runs

![MN5 LFP proxy: before L5 tuning vs Option 2](out/mn5_before_vs_option2_comparison.png)

Both runs are 5031 cells and 200 s on MN5; the plot shows 20–30 s, after the
start-up period, on identical y-scales. For each run the upper trace is the
0.5–2 Hz band plus the 10–15 Hz band, the lower trace is the 10–15 Hz band
alone, and grey boxes are detected spindle events (question 1).

- **Top, job 45171023:** before the L5 tuning. It predates three changes: L5
  recurrence, weaker L5 inhibition (`g_i_e_l5`) and slower L5 adaptation
  (`taur_l5e_rs`).
- **Bottom, job 45453403:** Option 2, with all three changes (L5 NMDA
  recurrence, `tau2 = 200 ms`, `taur = 120 ms`, `g_i_e_l5 = 0.02`).

Option 1 has not been run on MN5, so this is not the Option 1 vs Option 2
comparison from question 1, and the difference reflects the three L5 changes
together, not the recurrence alone.

Over the full run (first second excluded):

| | Spindle events | Rate | Interval between events | CV of intervals | Median event length | SO amplitude (RMS) |
|---|---:|---:|---:|---:|---:|---:|
| Before L5 tuning | 279 | 1.40/s | 712 ± 2 ms | 0.00 | 106 ms | 0.090 |
| Option 2 | 258 | 1.30/s | 771 ± 384 ms | 0.50 | 85 ms | 0.119 |

Before the tuning the network is a metronome: one spindle every 712 ms with
±2 ms jitter and a steady sawtooth slow wave. Option 2 has irregular spindle
timing, with bunched events and long gaps, and a larger slow wave with deep
troughs (around 21 s and 26 s). The "grey box = amplitude peak, not a full
spindle" caveat from question 1 applies to both.

Regenerate with:

```bash
python3 neuron/ctx_analyze.py --compare "MN5 job 45171023: before L5 tuning (no L5 recurrence)=res/2026-08-31/ctx_nrn_45171023.npz,MN5 job 45453403: Option 2 (L5 NMDA recurrence; tau2 200 ms / taur 120 ms)=res/2026-09-07/ctx_nrn_45453403.npz" --outdir out --tag mn5_before_vs_option2 --window-start 20000 --window-len 10000
```

Labels can't contain `,` or `=`, because `--compare` splits on them.

### How this differs from a real LFP

- **A layer's trace is built from that layer's own output spikes.** A real LFP
  comes from the synaptic currents flowing into the cells at the recording
  site. The proxy is really "the current these spikes would cause downstream",
  credited to the layer that fired them.
- **One 100 ms kernel is used for every spike,** but the model's AMPA synapses
  decay in 2 ms and GABA_A in 6–8 ms. Only the optional L5 NMDA synapses are
  that slow (115–200 ms). So the proxy is much smoother than the real synaptic
  currents and can add slow-wave content that doesn't exist in them. The
  docstring explains the history: the kernel used to be 10 ms / 60 ms, which
  hid the slow content. Moving to 100 ms fixed that, but may now overstate it.
- **Normalising each layer gives all four equal weight,** whatever their cell
  counts and firing rates.
- **The model has no electrode position and no dendritic geometry.** The
  deep-layer sign flip is a convention, not a calculated dipole.
- **Inhibition contributes nothing,** although real LFPs carry a large
  GABAergic component, especially during spindles.

### More standard options, if needed

- **Record synaptic currents in NEURON** (the `Exp2Syn`/`NMDA` currents per
  population) and use a current-based proxy. Mazzoni et al. (2015, *PLoS
  Comput Biol*) showed that a weighted sum of AMPA and GABA currents tracks the
  real LFP well.
- **Calculate the extracellular potential directly** with NEURON's
  `extracellular` mechanism or LFPy. Both need cell morphology, which the
  current point-like cells mostly lack.

## 4) Can we calculate the EEG from the LFP? If yes, how?

Not from the LFP trace we have now. But yes from the same simulation, by
computing a different quantity.

### Why the current LFP trace can't be converted

An LFP and an EEG are two measurements of the same thing: currents crossing
neuron membranes. They measure it in different ways.

- **An LFP** is the voltage at one point inside the tissue, dominated by nearby
  cells.
- **An EEG** is the voltage on the scalp. At that distance, the tissue's
  currents only matter through their net **current dipole moment**: roughly,
  how much current flows along the pyramidal cells' long axis. The scalp signal
  is that dipole seen through the brain, cerebrospinal fluid, skull and skin.

The standard chain is therefore:

> membrane currents → current dipole moment **p(t)** → head model (volume
> conductor) → scalp voltage

This uses the dipole, not the LFP. Our current LFP proxy (question 3) doesn't
contain the dipole: it's built from spike times with one standard kernel and a
sign convention we chose, not from real currents.

One useful consequence: the head model is quasi-static, meaning it has no
delays or frequency filtering. For one cortical patch with a fixed
orientation, the scalp signal is just the dipole's time course multiplied by a
constant. So **the EEG waveform is the dipole waveform**. The head model only
sets the amplitude in µV and how the signal spreads across electrodes.

### How to do it in this project (best option first)

**A. Compute the dipole directly in NEURON (most accurate, and practical
here).** Our pyramidal cells (`PYCell` in `neuron/cortex_neuron.py`) have a
25 µm soma plus a 300 µm dendrite split into 5 segments, and synapses sit at
the middle of the dendrite. That is enough geometry to compute a dipole:

1. Turn on NEURON's fast membrane-current output (`cvode.use_fast_imem(1)`).
2. For each cell, compute `p_z(t) = Σ_seg i_membrane(seg) · z(seg)`, with z
   measured along the soma–dendrite axis. Assume all pyramidal dendrites are
   parallel and point towards the brain surface.
3. Sum over the cells on each MPI rank, add up the ranks, and save at 1 kHz
   next to the spikes.
4. Which cells to include:
   - Include L2/3E, L5E and L6E.
   - Leave out L4E, which in reality is mostly stellate cells whose currents
     cancel out.
   - FS interneurons have no dendrite, so they contribute nothing anyway.
   - Leave out the thalamus: it is deep and its geometry cancels out, so it
     isn't visible in EEG.

This needs code changes in `neuron/ctx_thalamus_mpi.py` and new MN5 runs.

**B. Current-based EEG proxy.** Martínez-Cañada et al. (2021, *PLoS Comput
Biol*) built an EEG proxy for simple point-neuron networks. It is a weighted
sum of the AMPA and GABA currents onto pyramidal cells, and was validated
against full biophysical models. It needs synaptic currents to be recorded, so
it also needs new runs. It's simpler than A but less faithful to our actual
geometry.

**C. From spikes only, using the existing `.npz` files.** The kernel method of
Hagen et al. (2022, *PLoS Comput Biol*) convolves each population's firing rate
with a kernel for each projection, to predict the dipole. Our current proxy is
a crude version of this, with one kernel for everything. Using per-projection
kernels with the right synapse time constants and signs would already help.
This is the least accurate option, but it needs no new simulations.

**Then the head model**, only needed for µV values or electrode maps: LFPy has
a four-sphere head model (Næss et al. 2017) and the realistic New York Head
model (Huang et al. 2016). Hagen et al. (2018) describe using them in LFPy 2.0.

### Caveats

- **Amplitude.** 5031 cells is a tiny patch. A scalp EEG needs several cm² of
  cortex active together, which is millions of neurons. Absolute µV values
  therefore depend on assuming how many columns act together. The total scales
  roughly with the number of columns if they are perfectly in sync, or with its
  square root if they are independent. The waveform, spectrum and
  spindle/slow-oscillation timing are what can be compared with real data, not
  the amplitude.
- **Geometry is simplified.** Every synapse sits at the same point on a single
  straight dendrite. So the dipole's sign mostly reflects excitation versus
  inhibition at that spot. Real L5 cells, for example, also get input on their
  apical tuft.
- **The benefit.** With an EEG-like signal, we could apply standard sleep-EEG
  spindle detectors and compare with published values: spindle density per
  minute at central electrodes, and how spindles line up with slow
  oscillations. That's a better target than the LFP-proxy comparisons so far.

### References

- Hagen E, Næss S, Ness TV, Einevoll GT (2018). Multimodal modeling of neural
  network activity: computing LFP, ECoG, EEG, and MEG signals with LFPy 2.0.
  *Front Neuroinform*.
- Hagen E et al. (2022). Brain signal predictions from multi-scale networks
  using a linearized framework. *PLoS Comput Biol*.
- Huang Y, Parra LC, Haufe S (2016). The New York Head — a precise
  standardized volume conductor model for EEG source localization and tES
  targeting. *NeuroImage*.
- Martínez-Cañada P, Ness TV, Einevoll GT, Fellin T, Panzeri S (2021).
  Computation of the electroencephalogram (EEG) from network models of point
  neurons. *PLoS Comput Biol*.
- Næss S et al. (2017). Corrected four-sphere head model for EEG signals.
  *Front Hum Neurosci*.
