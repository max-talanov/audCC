# Plan: the stimulation experiments requested by Edo (30 Sept 2026)

Source: "from Edo 30 September 2026" (two pages: conditions to vary, measures, Graph A).
Builds on `PLAN-auditory-input.md` (C1/C2/C3 done) and the production NEURON model
(`neuron/ctx_thalamus_mpi.py`, 5,031 cells at scale 1.65, SO 0.77 Hz with nested spindles).

## 1. What Edo asks for

**Conditions (what we vary)**

| # | Parameter | Levels |
|---|---|---|
| P1 | SO phase at stimulation | rising, peak, falling, trough (4 bins of 90 deg) |
| P2 | Spindle at stimulation | absent / present |
| P3 | If present: spindle phase relative to SO | 4 bins (= "phase of SO at the spindle's peak") |
| P4 | Number of stimuli | 1 or 2 |
| P5 | Duration | to choose (proposal: 25, 50, 100, 200 ms) |
| P6 | Intensity | to choose (proposal: 4 levels, 30-75 dB) |

**Measures (stim vs no stim)**

- M1 Event-related potentials locked to the stimulus, per cortical layer and in the thalamus (Edo's question to Max: yes, see 3.1).
- M2 Event-related activity per neuron type, locked to the stimulus.
- M3 SO power after stimulation (0-1 s window).
- M4 P(spindle) after stimulation (0-1 s).
- M5 Spindle power after stimulation (0-1 s).
- M6 Phase of the spindle relative to SO after stimulation.
- **Graph A:** x = SO phase at stimulation, y = spindle-to-SO phase, colour = effect (M3-M6 each get one such graph).

## 2. What already exists and what is missing

| Need | State |
|---|---|
| Tones at chosen times, rank-independent, sham arm bit-identical to the plain run | Done (C1/C2, C3). Sham = the exact unperturbed counterfactual: the effect of a tone is stim minus sham *at the same time*, with the phase taken from the sham. This is the main tool of this plan |
| Closed loop with UP detector + delay | Done (C3, 3 seeds on MN5). Delay d is only a crude phase proxy (SO is irregular, CV ~0.4) |
| Per-level LFP, spindle detection (RE volleys), UP onsets, P(UP)/P(spindle)/TC rate vs d | Done (`ctx_analyze.make_levels_figure`, `volley_stats`, `aud_phase`, `cl_analyze`) |
| SO **phase** (continuous) and spindle **phase relative to SO** | Missing |
| Spindle **presence** at the tone time | Partly (`aud_phase` has UP/DOWN; spindle trains exist offline) |
| Two stimuli with an inter-stimulus interval | Missing (one `queue_tone` per detection) |
| Duration / intensity sweeps | `--cl-dur`, `--cl-level` exist, but no dose-response calibration |
| SO power, spindle power, spindle phase after stimulation, Graph A | Missing |
| Thalamic LFP | Missing (only spikes; a proxy must be defined) |

**Known blocker:** the thalamus is too excitable (evoked spindles 58-94% even at 30 dB in C1/C2, and
1.7 TC spikes/cell in 50 ms in C3). Intensity and duration sweeps are meaningless until a dose-response
exists, so WP1 comes first.

## 3. Work packages

### WP0 - Definitions to agree with Edo (before any run)

1. Phase convention: 0 deg = which point? Proposal: cortical composite LFP (L2/3+L4+L5+L6), SO band 0.3-1.5 Hz, Hilbert phase, with peak = maximum UP-state activity (depolarised/firing), trough = DOWN-state minimum; rising = DOWN to UP transition, falling = UP to DOWN. In the model's LFP proxy the sign must be fixed once and documented.
2. Spindle phase: SO phase at the spindle's envelope maximum (as written on page 1). For a spindle still in progress at the tone this peak may lie in the future, so classification is done after the fact from the sham counterfactual. State this in the paper methods.
3. Two-stimulus protocol: inter-stimulus interval (Ngo-style 1.075 s? other?). Which phase defines the pair: the first stimulus.
4. Analysis window for M3-M6: 0-1 s as written. A 1 s window holds about 0.8 SO cycles, so a spectral "SO power" is poorly estimated. Proposal: SO-band amplitude (RMS of the 0.3-1.5 Hz filtered LFP) over 0-1 s, plus 0-2 s as a robustness check. Edo to confirm.
5. "Hongi's paper" (ERP reference): need the reference to match the ERP layout (bins, baseline, polarity).
6. Duration and intensity levels (WP3).

### WP1 - Calibrate the thalamus (gate for P5/P6)

- Intensity (dB) to response curve: TC spikes/cell, evoked-spindle probability, evoked-UP probability, at 5 levels, sham-controlled, 2 seeds.
- Target: a graded response (evoked spindle probability from ~10% to ~60% across the range), not saturation.
- Levers, in order: lower the dB range; stronger K_Na refractoriness (3e-4, gain 0.6 was only validated in the isolated thalamus; run it at full scale first); IC to TC synaptic weight (`auditory_input.py`).
- Deliverable: chosen intensity levels and the dose-response figure. Cost: ~10 runs.

### WP2 - Analysis code (local, no cluster time)

New module `neuron/stim_phase.py` (offline, from the sham `.npz`):

1. `so_phase(npz)`: composite cortical LFP, band-pass, Hilbert, phase at every ms. Validate: phase distribution vs UP onsets (UP onset should fall at the "rising" bin).
2. `spindle_events(npz)`: onset, end, envelope peak from RE volleys (`volley_stats`); spindle count by spikes, not LFP, because the LFP spindle band partly reflects filter ringing.
3. `classify_tone(t)`: SO phase bin, spindle present/absent, spindle phase bin. Each arm is classified from **its own pre-tone data**, because the arms diverge after the first tone and the sham is no longer the counterfactual at later tones. This needs a causal phase estimate (an acausal filter would leak the tone's response); `stim_phase.py` has one and validates it against the acausal phase on sham runs.
4. Measures M3-M6 per tone and per arm (stim at t, sham at t): SO amplitude, P(spindle), spindle power (RE/TC rate in 10-15 Hz plus cortical LFP), spindle phase after the stimulus; the effect is the paired difference stim minus sham.
5. `stim_figures.py`: ERPs per layer and thalamus by SO phase bin (M1), PSTHs per cell type (M2), Graph A heat-map plus 3D-bar version for M3-M6, with counts per cell printed and cells with n below threshold greyed out.
6. Thalamic field potential: not a measured LFP. Define a proxy (TC+RE summed synaptic current, or filtered population rate) and label it as such.

Tests: unit tests on synthetic sinusoid+burst data; check on the existing 3 `c12_sham` and 3 `cl_*` runs that classification reproduces the known ~25% spindle occupancy and ~68% UP time.

### WP3 - Stimulus protocol in the simulator (done 2026-10-16; E0 prepared in `MN5_NEURON.md`)

Small changes in `closed_loop.py`, `auditory_input.py`, `ctx_thalamus_mpi.py`:

- `--cl-n {1,2}` and `--cl-isi <ms>`: second tone queued at t+isi (same deterministic fibre-train rule, event index distinct).
- Open-loop schedule mode `--rand-tones <min_gap>:<max_gap>`: tone times drawn from (seed, index), no detector, rank independent. This is the main delivery mode (see section 4).
- Duration and level already parameters; log them per tone.
- Keep the sham rule: sham = same schedule, no sound, must equal the plain run bit-for-bit (re-check on 100 ranks; it was verified for C3).

### WP4 - Runs on MN5

See section 5 for the matrix. Same recipe as C3 (`MN5_NEURON.md`): upload `neuron/` with `rsync --delete`, check `closed_loop.py` is present, submit stim/sham pairs, verify bit-identity of the sham arms and the ~1,190 s wall time.

### WP5 - Analysis, figures, report

Graph A for M3-M6, ERP and PSTH figures, effect tables with CIs, a short results section in `PLAN-auditory-input.md` and the README figure, a methods paragraph for Edo.

## 4. Delivery mode: open-loop random tones, phase assigned after the fact

Why not a phase-locked closed loop as the main mode:

- SO is irregular (CV ~0.4), so an online phase predictor has large error; calibrating one is a project of its own.
- Edo's Graph A needs the whole phase range *and* all spindle relationships, including ones that happen rarely. Random tones sample the natural distribution; phase and spindle state are then read from the sham.
- ~~The sham is exact, so no phase estimation noise enters the classification.~~ Corrected during WP2: only the first tone of a run has an exact sham counterfactual. Later tones are classified with a causal phase estimate: median error 25 deg, 64% of tones in the correct 90-deg bin (3 sham runs; trough tones are often called rising because the DOWN-to-UP transition is unpredictable). Matched sham tones are still useful as the control arm, classified the same way.

Constraints on the schedule:

- The spindle refractory period is 4-8 s (K_Na), and an evoked response lasts about 1 s. Use gaps drawn uniformly from 5-9 s, so ~57 tones per 400 s run.
- The first 2 s are skipped (as before).
- Arms diverge after the first tone, so stim tones are compared with sham fake tones of the same class (SO phase bin, spindle state), not tone by tone.

A phase-targeted closed loop (online SO phase predictor, plus spindle onset detector) is a later top-up for cells that stay empty. Its feasibility is evaluated offline first on the sham runs: circular error of predicted vs true phase.

**Geometric caveat to tell Edo:** a spindle lasts ~0.5 s, about 140 deg of an SO cycle, and in this model spindles sit in the UP state. The (SO phase at stimulation, spindle-to-SO phase) grid is therefore mostly populated along a diagonal band, and off-diagonal and off-UP cells will stay nearly empty with natural spindles. Only stimulus-evoked spindles can populate them, and that is itself a result. Graph A will show empty cells explicitly instead of interpolating.

## 5. Run matrix and budget

Per run: 400 s simulated, 100 ranks (about one node), ~1,190 s wall, ~0.33 node-hours. Each experiment needs a stim/sham pair, so 2 runs per seed.

| Block | Content | Seeds | Runs | Tones (stim) | Node-hours |
|---|---|---|---|---|---|
| E0 calibration | intensity dose-response (WP1) | 2 | ~10 | ~100 | ~3 |
| E1 reference | n=1, 50 ms, mid intensity, random phase | 16 | 32 | ~900 | ~11 |
| E2 one-factor sweeps | duration x3 extra, intensity x3 extra, n=2 with 2 ISIs; one factor at a time | 8 each | 8 combos x 16 = 128 | ~3,600 | ~43 |
| E3 top-up (only if cells stay empty) | targeted tones | tbd | tbd | tbd | tbd |

About 57 node-hours for E0-E2, to be checked against the MN5 allocation. E1 gives ~900 tones: with ~25% spindle occupancy and 4 SO bins that is ~55 per SO bin with a spindle present, so the full spindle-phase grid is underpowered at the single-cell level. Options if more power is needed: coarser bins (2 x 2), pooling the neighbouring cells, or more seeds (each extra 16-run block adds ~900 tones for ~11 node-hours).

## 6. Statistics

- Paired effect per tone (stim minus sham at the same time). Continuous measures: bootstrap over tones, resampling seeds as blocks. Binary (spindle yes/no): Wilson CIs and Fisher, as in `cl_analyze`.
- Phase measures (M6): circular mean and resultant vector length per cell; Rayleigh test; Watson-Williams (or bootstrap) between stim and sham.
- Cells with fewer than 15 tones are marked, not shown as estimates.
- Multiple comparisons: Graph A cells are descriptive; report the omnibus test (phase bin x stim) first.

## 7. Order of work and milestones

1. WP0 questions answered by Edo (blocking for definitions only; WP2 can start with the proposals).
2. WP2 offline analysis on existing runs (about a week); first Graph A from the C3 data, which has 439 tones, as a dry run with coarse bins. Useful for the discussion with Edo even before new runs.
3. WP1 calibration run (1 submission).
4. WP3 code + bit-identity check on MN5, then E1.
5. E2 sweeps, final figures and report.

## 8. Risks

| Risk | Mitigation |
|---|---|
| Thalamus saturates, no dose-response | WP1 first; K_Na at full scale; lower dB range |
| Too few tones per phase x spindle cell | coarser bins, more seeds, E3 targeted tones |
| SO power in a 1 s window is noisy | amplitude measure, 0-2 s check |
| LFP "spindle band" overstated by filter ringing | spindle measures from RE/TC spikes, LFP only as a check |
| Thalamic LFP proxy not comparable to a measured field | label as proxy; show also the TC/RE rate |
| Evoked responses overlapping in time | gaps of 5-9 s; skip tones inside the previous response window |
| Model SO phase convention differs from the experiment | WP0.1, documented once in the methods |

## 9. Stimulation rates from the literature (added 2026-10-23)

Focus 2 of the current work: compare stimulation frequencies used in human
sleep studies, in the model. What the literature search gave (only what was
confirmed in the sources; the rest is marked "not verified"):

| Protocol | Source | Rate / structure |
|---|---|---|
| Rhythmic, fixed rate | Ngo et al. 2013 (J Sleep Res), [arXiv 1407.7999](https://ar5iv.labs.arxiv.org/html/1407.7999) | 0.8 Hz tones, compared with random tones and sham; entrains SO once sleep is established, delays sleep onset if started before sleep |
| Rhythmic, fixed rate, in blocks | Tononi group (Bellesi et al. 2014, as described in [Frontiers 2017](https://www.frontiersin.org/journals/human-neuroscience/articles/10.3389/fnhum.2017.00109/text)) | 50 ms tones at a fixed 1 Hz, 15-tone blocks alternating with silence |
| Phase-locked pair | Ngo et al. 2013 (Neuron) | clicks in phase with the SO up state enhance SO, spindles and memory; out of phase ineffective. The 1.075 s click interval is **not verified** from the sources found |
| Phase-locked, more clicks | Ngo et al. 2015 (J Neurosci), as cited in [Debellemanière et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9788124) | little added benefit beyond two clicks, attributed to a refractory period, especially for fast spindles |
| Phase-locked 5-click sequence | Papalambros et al. 2017, same citation | a five-click sequence benefited SO and fast spindles for each click |
| Trains with SO-matched interval | Debellemanière et al. 2022 | 10 clicks of 50 ms, 40 dB; only the first is closed-loop, the others follow at a fixed interval matched to each subject's mean SO frequency; 10 s pause after each train |
| Single phase-locked stimuli | Koller et al. 2025 ([bioRxiv](https://www.biorxiv.org/content/10.1101/2025.03.04.639851.full.pdf)) | one 286 ms sound per stimulation, timed to the SO negative-to-positive transition (0 deg in their convention); 4.5 s pause after each; SO and delta rise, **fast spindles fall**, partly because the sound interrupts ongoing spindles; the response to later pulses of a train "quickly fades out" |

Open point for Edo: "Hongi's paper" is probably by Hong-Viet Ngo (the Ngo
papers above) but that is a guess; which paper's ERP layout should we copy?
Koller et al. use ERP epochs from -1.5 to 3.5 s around the stimulus.

### What the model can compare

The model's SO rate is 0.77 Hz (period about 1.3 s), so 0.8 Hz stimulation is
almost resonant with it, 1 Hz drifts through all phases, 0.5 Hz is slower than
the SO and 1.67 Hz is about twice as fast. The questions: does a rhythmic train
entrain the UP states (vector strength of UP onsets in the stimulation period),
does the response fade along the train (the refractory period of the
literature), does the UP and spindle rate change during and after the train.

**F1 design** (`MN5_NEURON.md`, "F1 frequency comparison"): 50 ms stimuli at 35
dB, trains of n = 1 (the single-stimulus reference) or n = 10, stimulation rates
0.5 / 0.8 / 1.0 / 1.67 Hz (`--ol-isi 2000:1250:1000:600`), 5-9 s after the end
of each train. 8 seeds x (stim, sham) = 16 jobs, about 32 trains per rate. The
level follows E0: TC response is graded (0.62 to 2.38 spikes/cell from 25 to 75
dB), the evoked-spindle probability is already 60-80% at 25-75 dB, so 35 dB is
low but not subthreshold; the K_Na variant is a later control.

Analysis: `neuron/stim_train.py` (per pulse index: TC response, P(UP),
P(spindle); per train: UP rate, entrainment vector strength with Rayleigh test,
UP and spindle rate 1-6 s after the train; stim vs sham), tested on synthetic
data (`neuron/test_stim_train.py`).
