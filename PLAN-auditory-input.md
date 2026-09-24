# Plan: auditory input to the thalamocortical loop, awake and asleep

Status: **Stage A done (2026-09-24, local, see "Stage A results"); Stages
B and C not started.** Builds on the production
NEURON model `ParallelCorticoThalamicNet` (`neuron/ctx_thalamus_mpi.py`;
architecture in `Edu-questions.md` Q2). Runs alongside `PLAN-spindels.md`,
and the spindle-related experiments here depend on it (see "Dependencies").

## Goal

Add a sound-driven input to the auditory thalamocortical column and one
brain-state setting that switches the same network between:

- **Wake:** thalamus in tonic (relay) mode, cortex asynchronous, sounds
  relayed roughly linearly, MGB → L4 → L2/3.
- **NREM sleep:** the current regime. Hyperpolarised thalamus in burst mode,
  a cortical slow oscillation (SO), and, once `PLAN-spindels.md` succeeds,
  spindles. Sounds arrive at different SO / spindle phases.

The scientific output is how the **same stimulus** is transmitted, gated or
turned into a sleep event (evoked K-complex / slow wave, evoked spindle)
depending on state and phase. It ends with the closed-loop acoustic
stimulation experiment (Ngo et al. 2013) that Mushtaq et al. 2024, the source
of this model's parameters, studies.

## What the model has now (relevant facts)

| Aspect | Now | Consequence for this plan |
|---|---|---|
| External input | **None.** No IC/brainstem source, no background noise (`PLAN-spindels.md` diagnosis 1) | An input population must be added from scratch |
| Sleep state | Set only by hyperpolarised `e_pas` in `_make_cell` (TC −80, RE −82, cortex −70 / FS −68 mV) | There is no knob to leave sleep; wake must be added as a new mechanism |
| Cortical SO | Intrinsic: L5 `PYCellIB` (I_NaP + slow SK2) + IB gap junctions, ~1.4 Hz, optional L5 NMDA recurrence (Options 1/2) | In wake these mechanisms must be turned down, not removed |
| TC → cortex | TC → L4E only (`_wire_thalamocortical`); no TC → L4I | No feedforward inhibition, so evoked responses won't be sharp. The NEST model has `thalamus_E → L4_I`; the NEURON model lacks it |
| Topography | None. Fixed-convergence random draws (`_draw`) | No tonotopy. `PLAN-spindels.md` Stage 1 step 2 adds topographic RE ↔ TC wiring, so both plans should share one coordinate |
| Synapses | `Exp2Syn`, no short-term plasticity | No adaptation to repeated sounds (SSA) without adding depression |
| Spindles | Single RE volleys, not spindles (`Edu-questions.md` Q5) | Spindle-gating experiments must wait for `PLAN-spindels.md` Stage 2 |
| Determinism | Everything is a function of gid, independent of rank count | Stimulus spike trains must follow the same rule |

## Design

### D1. The input population: lemniscal IC → MGBv drivers

- **Source:** artificial "IC fibres" (inferior colliculus, central nucleus),
  not cells. Each fibre has a spike train computed from the stimulus,
  **seeded by fibre id**, so every rank can compute any fibre's train without
  MPI traffic. This keeps results identical at 2, 4 or 100 ranks.
- **Delivery:** `NetCon(None, syn)` onto a dedicated AMPA synapse per TC cell;
  spikes are queued with `nc.event(t)` from an `FInitializeHandler` (open
  loop) or at run time (closed loop, D5). One mechanism serves both modes.
- **Driver synapse, not modulator:** each TC gets few, strong inputs (k ≈ 5–10
  fibres, unlike the k = 100 weak L6 inputs), and optionally a slow NMDA part
  (`nmda.mod` already exists). Calibrate so that 1–2 coincident fibres fire a
  **depolarised** (wake) TC cell. Whether they fire a **hyperpolarised** (NREM)
  cell is what the experiments measure; don't tune for it.
- **Rate model per fibre:** `r(t) = r_spont + A(level) · tuning(f_stim, CF) ·
  [onset transient + sustained part]`, inhomogeneous Poisson with a 1 ms
  refractory period. Onset ~10–20 ms, then decay to a lower sustained rate.
  Keep every parameter explicit in one stimulus spec.
- **Stimulus types:** pure tone (frequency, level, duration), click, click
  train (e.g. 40 Hz, for the auditory steady-state response), noise burst,
  oddball sequence (standard/deviant) for stimulus-specific adaptation (SSA).
- **Optional later:** GABAergic IC → MGB (tectothalamic) inhibition.

### D2. Tonotopy, shared with the spindle plan

- Give TC, RE and L4 (then L2/3, L5, L6) a position `x ∈ [0, 1]` on one
  tonotopic axis. Channel = bin of `x`. The same `x` is the ring/line
  `PLAN-spindels.md` Stage 1 step 2 needs for local RE ↔ TC wiring, so
  **implement it once, in `_draw`, behind a flag**.
- Stage A starts with **one channel** (the whole column shares one
  characteristic frequency, CF). The stimulus then only has level, timing and
  type. Tonotopy (8–16 channels, Gaussian tuning on IC → TC, TC → L4, and
  broader L4 → L2/3) is Stage B2.
- Honest scope: one column is a narrow slice of A1. A tonotopic model is a
  **strip** of A1 made from the same cells, not a single column any more.

### D3. Thalamocortical transmission

- Add **TC → L4I** (feedforward inhibition), as in the NEST model and
  Mushtaq 2024 Table 3 (TC → IN). Without it evoked L4 responses have no
  sharp window.
- Add **short-term depression** (Tsodyks–Markram; a `tmgsyn.mod` is available
  on ModelDB) on IC → TC and TC → L4, behind a flag. It is needed for SSA and
  for state-dependent transmission: the tonic spontaneous rate in wake keeps
  TC → L4 synapses depressed, while bursts after silence in NREM hit rested
  synapses (Swadlow & Gusev 2001).

### D4. Brain state as neuromodulation, not new wiring

The switch from NREM to wake is modelled as the neuromodulatory action of
ACh / NE / 5-HT / histamine (McCormick 1992; Bazhenov et al. 2002). One state
object sets all of it; the wiring never changes.

- **K⁺ leak (`kleak.mod`, new):** `i = g_kl · (v − E_K)`, `E_K ≈ −95 mV`, on
  every cell. It is the Hill & Tononi / Bazhenov 2002 sleep-depth knob the NEST
  model already uses (`g_KL`). Neuromodulation closes it in wake, which
  depolarises TC out of burst mode (I_T inactivated) and depolarises cortex.
  Knobs `g_kl_tc`, `g_kl_re`, `g_kl_cx`, separately, since ACh depolarises TC
  but hyperpolarises RE while NE depolarises both.
- **Recalibration step (first):** split today's `e_pas` into `pas + kleak`
  such that **NREM rest potential and input resistance are unchanged**
  (combined reversal `(g_pas·e_pas + g_kl·E_K)/(g_pas + g_kl)` = today's
  value). NREM must then reproduce the current spike counts (regression gate,
  Stage A).
- **Cortical adaptation:** ACh blocks the slow AHP. Scale `gkbar_sk2` in
  cortex (and `PYCellIB`'s slow `taur`) by a wake factor (~0.3, to sweep).
  This is what removes the SO in wake.
- **L5 IB drive:** scale `gnap` on `PYCellIB` and optionally the IB gap
  conductance `g_l5_gap` (a model construct used for SO synchrony; see risk
  R2).
- **Synaptic modulation (optional, Stage B):** ACh reduces intracortical
  excitatory transmission and strengthens thalamocortical transmission
  (Gil et al. 1997; Hasselmo 1995): scale `g_ff` / `g_l5_rec` down and
  `g_tc_l4` up.
- **Background synaptic noise (new):** per-cell excitatory + inhibitory
  Poisson input (`NetStim`, `noise = 1`, Random123 streams keyed by gid).
  Wake needs it for an asynchronous, high-conductance cortex (Destexhe et al.
  2003). In NREM keep it **off by default** until `PLAN-spindels.md` decides
  whether noise is also its desynchronising lever (its diagnosis 1).
- **Time-varying state (Stage C3):** every state parameter can follow a
  schedule (`Vector.play` on the RANGE variables). That gives wake → NREM
  transitions, and a slow ~0.02 Hz infraslow modulation (review §3.2 /
  `docs/spindle_review_mapping.md` §5.1B) if wanted.

State presets live in one file, e.g. `neuron/brain_state.py`:

```python
STATES = {
    "nrem": dict(g_kl_tc=..., g_kl_re=..., g_kl_cx=..., sk2_scale=1.0,
                 gnap_scale=1.0, l5_gap_scale=1.0, bg_rate_e=0.0, bg_rate_i=0.0),
    "wake": dict(g_kl_tc=0.0, g_kl_re=..., g_kl_cx=..., sk2_scale=0.3,
                 gnap_scale=..., l5_gap_scale=..., bg_rate_e=..., bg_rate_i=...),
}
```

The `nrem` values come from the recalibration; the `wake` values from Stage B1.

### D5. Closed-loop stimulation (online SO phase)

- Run `psolve` in chunks (e.g. 5 ms). After each chunk, `allreduce` the L5E +
  L6E spike count, feed a causal online SO estimate on rank 0 and broadcast a
  "stimulate at t" decision. Each rank then queues the fibre events for its own
  TC cells (D1). About 40 000 allreduces per 200 s run, which is cheap next to
  the simulation.
- Trigger rule (Ngo et al. 2013): detect a DOWN state (population silence
  after an UP), stimulate after a delay `d`. Sweep `d` to move the tone
  across SO phases.
- **Sham condition is mandatory:** the same detector, detections logged, no
  sound. Compare stimulated and sham at matched detection times.

### D6. Measurement

A new `neuron/aud_analyze.py`, reading the run's `.npz` plus the stimulus log
(which the run must save: tone onset times, type, level, channel, and
detector events).

- **PSTH** per population (and per channel), onset latency, evoked spike
  count in 0–50 ms and 50–300 ms windows.
- **TC firing mode:** burst vs tonic spikes, Lu et al. 1992 criterion
  (≥ 100 ms silence, then ≥ 2 spikes with ISI ≤ 4 ms). This is the direct
  readout of the state switch.
- **Cortical state:** ISI CV, pairwise spike-count correlation (50 ms bins),
  SO-band (0.5–2 Hz) power of the LFP proxy (`_synaptic_lfp_local`), fraction
  of time in population silences ≥ 100 ms (DOWN states).
- **Phase binning:** SO phase at each tone from the Hilbert phase of the
  0.5–2 Hz LFP proxy (offline, for open loop). Spindle on/off from
  `PLAN-spindels.md` Stage 0 (spike-based RE-volley trains), **not** a
  band-pass detector.
- **Evoked sleep events:** probability that a tone is followed within
  ~0.5–1 s by an SO cycle (evoked K-complex / slow wave) or a spindle train,
  against sham/no-tone windows at matched phase.
- **Tuning (Stage B2):** frequency response areas per layer and state.

## Dependencies on PLAN-spindels.md

| This plan | Needs from the spindle plan |
|---|---|
| Stage A, B (input, wake) | Nothing. Can start now |
| Stage C1, C2 (NREM input, SO phase, evoked slow waves) | Nothing. The SO exists today |
| Stage C2 spindle part, C3 closed loop on spindles | Stage 2 done (real spindles with the cortex attached) |
| D2 tonotopy | Shares the topographic coordinate with Stage 1 step 2; build it once |
| NREM background noise | Its decision on noise in NREM |

Every new mechanism is behind a flag with the old behaviour as default, so
the two plans don't break each other's runs.

## Stage A — infrastructure and NREM regression (local) — **done 2026-09-24**

- [x] `mod/kleak.mod`; add to every cell type; `nrnivmodl` build and the MN5
      upload list (`MN5_NEURON.md`).
- [x] Recalibrate NREM `e_pas` → `pas + kleak` (D4). **Regression gate:**
      `--state nrem` with no stimulus reproduces baseline spike counts per
      population within a few percent at `--scale 0.1`, same seeds.
- [x] `neuron/auditory_input.py`: stimulus spec (JSON, so MN5 needs no YAML;
      or a dict), fibre spike trains keyed by fibre id, IC → TC driver
      synapses, event queuing. One channel.
- [x] Unit checks (`neuron/aud_checks.py`): fibre trains identical at 1 and 4
      ranks; PSTH of the fibres alone matches the rate model; a single TC
      cell held at −68 mV vs −80 mV answers one fibre volley with tonic
      spikes vs an I_T burst. (−68, not −60: this `TCCell` fires tonically on
      its own from about −67 mV.)
- [x] `--state` and `--stim` flags on `ctx_thalamus_mpi.py` (and `STATE` /
      `STIM` in `run_ctx_nrn.sh`). No separate `--save-stim-log`: with
      `--stim` the schedule and all fibre spikes always go into the `.npz`.
- [x] `neuron/aud_analyze.py`: PSTH, TC burst/tonic classifier, cortical-state
      metrics (D6), plus `--compare` for regression and rank-count checks.

**Done when:** the NREM regression gate passes, and tones at `--scale 0.1`
produce a measurable TC response in NREM with the analysis script reporting
it.

### Stage A results

All at `--scale 0.1` (305 cells), 4 ranks, Option 2 flags. Raw output:
`res/2026-09-24/stageA_summary.txt`; figure:
`out/stageA_nrem60_tones_psth.png`.

- **Leak split is exact per cell:** legacy vs `nrem` voltage traces differ by
  < 1e-8 mV, same spikes, for TC, RE, L5 IB, RS pyramidal and FS cells.
- **Regression gate passes.** In the network, legacy and `--state nrem` give
  identical rasters for the first 4.9 s, then diverge (the network is
  chaotic, so a 1e-8 mV rounding difference grows). After 20 s the largest
  per-population count difference is 5.1% (L2/3E; TC 0.1%). The noise floor
  is larger: nudging every cell's `e_pas` by 1e-9 mV in the legacy model
  changes counts by up to 9.0% (L6E), divergence at 3.3 s. So the split
  changes nothing beyond the network's own sensitivity.
- **Rank count doesn't matter:** NREM + tones at 1 and 4 ranks give identical
  spike rasters.
- **Driver calibration:** w = 0.002 µS per fibre, a ~4 mV EPSP. At −68 mV two
  coincident fibres fire a TC cell and one doesn't; at −80 mV two fibres
  trigger a 4-spike I_T burst.
- **Tones reach the thalamus and L4 in NREM** (60 s, 20 tones of 60 dB,
  50 ms, vs a sham run without `--stim`; spikes/cell, Wilcoxon paired):

  | 0–50 ms | tone | sham | p |
  |---|---|---|---|
  | TC | 4.96 | 0.71 | 9e-5 |
  | RE | 2.10 | 1.00 | 0.07 |
  | L4E | 1.89 | 0.65 | 3e-4 |
  | L4I | 4.16 | 1.71 | 6e-4 |
  | L2/3E | 0.24 | 0.11 | 0.07 |
  | L2/3I | 2.08 | 0.78 | 0.015 |
  | L5E / L6E | 0.45 / 0.41 | 0.37 / 0.24 | n.s. |

  TC responds on 20/20 tones, median first-spike latency 7 ms after onset.
  Nothing changes in 50–300 ms, so no evoked slow-wave cycle is visible at
  this sample size (a Stage C question).

**Found along the way, for Stage B:**

- The placeholder `wake` preset (K⁺ leak fully closed, `E_OPEN["tc"] = −62`)
  turns an isolated TC cell into a **35 Hz pacemaker with no input**. Stage
  B1 has to calibrate it, starting from the fact that this `TCCell` fires
  tonically from about −67 mV.
- The current NREM network is less "bursty" than the plan's target table
  assumes: only **~25% of TC spikes are Lu-criterion bursts** (TC fires at
  ~13 Hz/cell at this scale), and pooled cortical E activity never goes
  silent for ≥ 100 ms, so the silence metric reads 0 in NREM. The
  wake/NREM contrast in B1 should be judged relative to this baseline, and
  the silence metric may need a per-layer or thresholded version.

## Stage B — the awake state (local, then MN5)

### B1. Find the wake operating point (single channel)

- [ ] Add TC → L4I and background noise (D3, D4).
- [ ] Sweep the wake knobs one at a time, in order: `g_kl_tc`/`g_kl_re`,
      `sk2_scale`, background rates, `gnap_scale`, `l5_gap_scale`. Keep each
      that helps, as in `PLAN-spindels.md`.
- [ ] Targets, **spontaneous activity**:

  | Measure | Wake target | NREM (for contrast) |
  |---|---|---|
  | TC burst fraction (Lu criterion) | low (< ~10 %) | high |
  | TC / RE mode | tonic single spikes | rebound bursts |
  | Cortical ISI CV | ~1 (irregular) | periodic, SO-locked |
  | Pairwise correlation (50 ms bins) | near 0 | high |
  | SO-band power (0.5–2 Hz) | ≥ 10× lower than NREM | peak near SO |
  | Population silences ≥ 100 ms | rare | every SO cycle |
  | Rates | low, sparse (E a few Hz, FS higher) | as now |

- [ ] Targets, **evoked (tones)**: short-latency, monotonic-in-level MGB
      response; L4 response a few ms after MGB; L2/3 after L4; an onset
      response larger than the sustained one; no evoked SO cycle. Check the
      order and relative latencies against rodent / primate A1 data rather
      than fitting absolute values.

**Done when:** both tables hold on ≥ 4 of 5 seeds at `--scale 0.1`, and NREM
still passes the Stage A regression gate with the new mechanisms at their
NREM values.

### B2. Awake auditory physiology

- [ ] Tonotopy (D2), 8–16 channels. Frequency response areas: narrow in MGBv
      and L4, broader in L2/3.
- [ ] Short-term depression (D3). Oddball sequence: deviant response larger
      than standard (SSA), stronger in cortex than in MGBv.
- [ ] 40 Hz click train: steady-state response in L4 and the LFP proxy.

**Done when:** tuning, SSA and 40 Hz following are present in wake. These are
qualitative sanity checks, not fits.

### B3. Full-scale wake run on MN5

- [ ] One 200 s wake run at `--scale 1.65` with a tone protocol; confirm no
      drift (window comparison, as in `Edu-questions.md` Q8).

## Stage C — auditory input during NREM

### C1. Open loop, SO phase (can start after Stage A)

- [ ] Tones at random times (inter-tone interval 2–5 s, jittered, so they
      sample all SO phases), several levels. 200 s per seed, several seeds.
- [ ] Measure, binned by SO phase (UP, DOWN, DOWN→UP, UP→DOWN):
  - TC transmission: burst vs tonic response, evoked spike count;
  - L4 and L2/3 evoked response vs the same tone in wake;
  - probability of an evoked SO cycle / K-complex; latency;
  - whether an evoked UP state drives the L6 → RE kick.
- [ ] Questions to answer, **measured, not assumed**:
  1. Is thalamic transmission reduced in NREM, and how much of it is a phase
     effect (DOWN vs UP)?
  2. Do tones in DOWN states evoke TC bursts that start an UP state
     (burst-mode "wake-up call", Sherman 2001)?
  3. Are early cortical responses (L4) preserved while later / L2/3 ones are
     reduced? Single-unit data suggest primary auditory cortex responses are
     largely preserved in NREM (Issa & Wang 2008; Nir et al. 2015), so
     "the thalamus closes the gate" is a hypothesis to test here, not a
     target.

### C2. Spindle phase (after PLAN-spindels.md Stage 2)

- [ ] Same protocol; additional bins: during a spindle, just after a spindle
      (refractory), no spindle.
- [ ] Measure: transmission during spindles vs outside (reduced in humans:
      Dang-Vu et al. 2011; Schabus et al. 2012); probability of a tone evoking
      a spindle; refractoriness after a spindle.

### C3. Closed loop and state transitions

- [ ] Closed loop (D5): stimulate at a delay after DOWN-state detection; sweep
      `d` across phases; sham control. Measure SO amplitude/regularity and
      spindle rate vs sham. Expected from Ngo 2013 / Mushtaq 2024: in-phase
      (UP) stimulation increases SO and spindle activity; out-of-phase does
      not. The spindle part needs C2.
- [ ] Wake → NREM transition: ramp the state over ~20–30 s with tones
      throughout; response vs time as TC switches from tonic to burst mode.
- [ ] Optional: tone level high enough to shift NREM toward wake
      (arousal) needs a stimulus → neuromodulator pathway the model doesn't
      have; list it as out of scope unless wanted.

**Done when:** C1 has a phase-by-phase table from ≥ 5 seeds at `--scale 0.1`
and one MN5 run at `--scale 1.65`; C2/C3 after the spindle plan allows.

## Out of scope for now

- The cochlea and brainstem: IC fibre rates come from the rate model in D1,
  not from a peripheral model. A cochlear front end (e.g. Zilany et al. 2014)
  could replace it later with the same fibre interface.
- Non-lemniscal MGB (dorsal/medial) and the matrix pathway to L1
  (`PLAN-spindels.md` also excludes it).
- Top-down / higher-order auditory areas, and an arousal pathway from sound
  to neuromodulator nuclei.
- REM sleep (a third state: aminergic off, cholinergic on). The same state
  object could express it later.
- The NEST model. It already has a light `auditory_rate` Poisson drive onto
  MGB; mirroring this plan there is not planned.

## Risks

- **R1. Recalibrating the leak shifts NREM.** Replacing `e_pas` with
  `pas + kleak` at equal rest potential and input resistance should be
  neutral, but the Stage A regression gate must pass before anything else.
- **R2. The SO machinery may resist wake.** The L5 IB population is held
  synchronous by gap junctions and slow SK2 (`neuron/README.md`, MN5 rounds
  1–3). If lowering SK2, `gnap` and background input doesn't give an
  asynchronous cortex, reducing `g_l5_gap` in wake is justified: those gap
  junctions are a model device for SO synchrony, not a documented cortical
  mechanism.
- **R3. Driver strength sets the result.** How strongly IC drives TC decides
  how much gets through in NREM. Calibrate it **in wake only** (D1) and keep
  it fixed across states, so the NREM results are predictions.
- **R4. Spindle-related results are premature** until `PLAN-spindels.md`
  Stage 2 passes. Before that, report only SO-phase effects.
- **R5. Trial counts.** Phase-binned statistics need many tones per bin:
  plan ~60 tones per 200 s run and several seeds per condition. At ~2×
  realtime on MN5 (100 ranks, 5031 cells) that is cheap.
- **R6. One column vs a tonotopic strip.** With tonotopy, each channel has
  fewer cells than today's whole column; check that per-channel SO and
  thalamic dynamics survive the split (the scale-dependence lessons in
  `neuron/README.md` apply).

## Files (planned)

| File | Change |
|---|---|
| `neuron/mod/kleak.mod` | new: K⁺ leak, the state knob |
| `neuron/mod/tmgsyn.mod` | new: Tsodyks–Markram depressing synapse (Stage B2) |
| `neuron/brain_state.py` | new: `nrem` / `wake` presets and schedules |
| `neuron/auditory_input.py` | new: stimulus spec, IC fibre trains, driver synapses, closed-loop hook |
| `neuron/ctx_thalamus_mpi.py` | flags `--state`, `--stim`, `--closed-loop`, `--tonotopy`; TC → L4I; background noise; stimulus log in the `.npz` |
| `neuron/aud_analyze.py` | new: PSTH, burst/tonic, cortical state, phase binning, evoked events |
| `MN5_NEURON.md`, `run_ctx_nrn.sh` | new files and flags |
| `neuron/README.md`, `Edu-questions.md` | results as they come |

## References

- Bazhenov M, Timofeev I, Steriade M, Sejnowski TJ (2002). Model of
  thalamocortical slow-wave sleep oscillations and transitions to activated
  states. *J Neurosci*.
- Dang-Vu TT, Bonjean M, Schabus M, et al. (2011). Interplay between
  spontaneous and induced brain activity during human non-rapid eye movement
  sleep. *PNAS*.
- Destexhe A, Rudolph M, Paré D (2003). The high-conductance state of
  neocortical neurons in vivo. *Nat Rev Neurosci*.
- Fernandez LMJ, Lüthi A (2020). Sleep spindles: mechanisms and functions.
  *Physiol Rev*.
- Gil Z, Connors BW, Amitai Y (1997). Differential regulation of neocortical
  synapses by neuromodulators and activity. *Neuron*.
- Hasselmo ME (1995). Neuromodulation and cortical function. *Behav Brain Res*.
- Hill S, Tononi G (2005). Modeling sleep and wakefulness in the
  thalamocortical system. *J Neurophysiol*.
- Issa EB, Wang X (2008). Sensory responses during sleep in primate primary
  and secondary auditory cortex. *J Neurosci*.
- Lu SM, Guido W, Sherman SM (1992). Effects of membrane voltage on receptive
  field properties of lateral geniculate neurons in the cat. *J Neurophysiol*.
- McCormick DA (1992). Neurotransmitter actions in the thalamus and cerebral
  cortex and their role in neuromodulation of thalamocortical activity.
  *Prog Neurobiol*.
- Mushtaq M, Marshall L, ul Haq R, Martinez T (2024). Possible mechanisms to
  improve sleep spindles via closed loop stimulation during slow wave sleep:
  a computational study. *PLOS ONE*.
- Ngo HVV, Martinetz T, Born J, Mölle M (2013). Auditory closed-loop
  stimulation of the sleep slow oscillation enhances memory. *Neuron*.
- Nir Y, Vyazovskiy VV, Cirelli C, Banks MI, Tononi G (2015). Auditory
  responses and stimulus-specific adaptation in rat auditory cortex are
  preserved across NREM and REM sleep. *Cereb Cortex*.
- Schabus M, Dang-Vu TT, Heib DPJ, et al. (2012). The fate of incoming
  stimuli during NREM sleep is determined by spindles and the phase of the
  slow oscillation. *Front Neurol*.
- Sherman SM (2001). Tonic and burst firing: dual modes of thalamocortical
  relay. *Trends Neurosci*.
- Swadlow HA, Gusev AG (2001). The impact of "bursting" thalamic impulses at a
  neocortical synapse. *Nat Neurosci*.
- Tsodyks MV, Markram H (1997). The neural code between neocortical pyramidal
  neurons depends on neurotransmitter release probability. *PNAS*.
- Zilany MSA, Bruce IC, Carney LH (2014). Updated parameters and expanded
  simulation options for a model of the auditory periphery. *J Acoust Soc Am*.
