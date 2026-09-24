# Plan: make the thalamus produce real sleep spindles

Status: **Stage 0 done; Stage 1 started (lever 1 tested, negative) —
2026-09-24.** See "Results so far" at the end of Stage 1. Background and evidence:
`Edu-questions.md` question 5, and the history in `neuron/README.md`.

## Goal and success criteria

A real spindle, as used in this plan:

| Property | Target | Now (MN5, 5031 cells) |
|---|---|---|
| Cycles per event (RE volleys 60–150 ms apart) | ≥ 6, typically 7–15 | 1 (second cycle at ~9% of the first) |
| Frequency within an event | 10–15 Hz | none (single volley) |
| Duration | 0.5–2 s | a few tens of ms |
| Envelope | grows, then fades (waxing/waning) | none |
| Recurrence | every few seconds, with a refractory period | every 0.6–0.7 s, locked to the cortical kick |
| TC participation | each TC cell fires every 2nd–3rd cycle; RE fires every cycle | every cell fires in every volley |

All criteria are measured **from spikes**, never from band-pass counts: a
10–15 Hz filter turns a single volley into ~3 cycles of ripple
(`Edu-questions.md` Q5).

## Diagnosis: why the loop fires only once

1. **Total synchrony.** All 91 RE cells fire within 50 µs (`neuron/README.md`,
   MN5 round 3); in the Option 2 run the median spread is 16 µs, with every
   cell in every volley. **Main cause, measured: the RE ↔ RE gap junctions
   are ~30× too strong** (`Edu-questions.md` Q7). Contributing: RE → TC is
   all-to-all (k = 91), every L6E cell fires ~6 ms before each RE volley, and
   the model has no noise. After one volley every TC cell is in the same
   refractory state, so no subset is left to carry cycle 2.
2. **No slow inhibition.** RE → TC is GABA_A only (8 ms decay) in
   `ctx_thalamus_mpi.py` (`_wire_re_tc`). GABA_B was tried in the serial model
   as a *linear* `Exp2Syn` (τ 60/200 ms), active on every spike: it acted as
   tonic hyperpolarisation and silenced TC. Real GABA_B needs a high-frequency
   RE burst (G-protein cooperativity), then gives a ~150–200 ms IPSP that
   re-primes TC's I_T.
3. **No waning or refractory mechanism.** Ca²⁺-dependent I_h (`ihca.mod`)
   exists but is off in production (`gh_tc = 0`).
4. **The cortex kicks too often.** Volleys arrive every 590–712 ms, which would
   cut off a 0.5–2 s spindle even if the thalamus could ring. Real slow
   oscillations are < 1 Hz, and only some UP states carry a spindle.

## Already ruled out (don't repeat)

From `neuron/README.md` (mostly the 10–40-cell serial model):

- loop gain (`g_tc_re`, `g_re_tc`) over a 7× range
- population size 10 → 40
- gap-junction strength alone — but only in the 10–40-cell model, where every
  cell is a near neighbour anyway; **not** ruled out at 91 RE cells, where it
  is the measured main cause of RE synchrony (Stage 1, step 1)
- heterogeneous resting potentials
- "progressive recruitment" — but only in a 10-cell network, where local
  wiring can't exist; **not** ruled out at 437 cells
- linear GABA_B (`Exp2Syn`) — **not** a test of real, cooperative GABA_B
- I_h at 20× the published conductance: moved the loop to 13 Hz, but events
  stayed < 0.25 s

The one lever with a measured positive effect: per-cell heterogeneity
(`het = 0.05`) tripled cycles per event (2.5 → 7.5) in the small model, but at
6–7 Hz. This points to phase dispersion as the missing ingredient.

## Stage 0 — a spike-based spindle measure

Before changing the model, make success measurable without the band-pass
filter.

- [x] Extend `neuron/volley_stats.py` with RE-volley **trains**: detect RE
      population volleys; consecutive volleys 60–150 ms apart belong to one
      train. (`--trains`; a volley = RE-rate peak with ≥ 10% of RE cells
      firing within ±20 ms.)
- [x] Report per train: cycles, duration, within-train frequency, and the
      fraction of TC cells firing per cycle; per run: trains per minute,
      inter-train interval. (Also RE first-spike SD per volley and a
      waxing/waning flag.)
- [x] Baseline both MN5 runs (45171023, 45453403). Expected: 1 cycle per train.
      **Result:** 279 and 357 trains, **100% single volleys**, RE and TC
      participation 100% per volley, RE first-spike SD 0.02 ms. The Stage A
      runs (305 cells, with or without tones): 98% single, max 2 cycles.
      `--self-test` recovers an 8-cycle 12 Hz waxing/waning train, a single
      volley and a 3-cycle 10 Hz train from a synthetic raster.
      Figure: `out/stage0_trains_vs_bandpass_20-30s.png` — every grey
      band-pass "spindle" coincides with a single-volley tick.

**Done when:** the measure reproduces "1 cycle per event" on both runs and
correctly counts cycles on a synthetic multi-cycle spike train.

## Stage 1 — isolated thalamus, local runs

Test the thalamus's own ability to ring, separate from cortical pacing. Each
run is TC + RE only, seconds of wall time, on a laptop.

- [x] New script `neuron/thal_ring_test.py`: TC + RE at MN5 sizes (346 / 91),
      same cell parameters as `ctx_thalamus_mpi.py` (`gsk_re = 1e-3`,
      `het = 0.05`, current synaptic values). **Don't reuse `tc_mpi.py`:** it
      still defaults to the old `gsk_re = 5e-5` that caused the RE runaway, and
      it uses a periodic drive.
- [x] Protocol: 1 s settle, then **one** brief excitatory kick to RE and TC
      (mimicking a single L6 volley), then ≥ 3 s free-running. Repeat for
      ≥ 5 seeds. (Also a **tone kick** through the auditory IC → TC input,
      `PLAN-auditory-input.md`; seeds vary cell heterogeneity via the new
      `het_seed`, the kick jitter and the tone.)
- [x] Baseline with current wiring. Expected: 1 cycle, reproducing MN5.
      **Result:** 1 cycle on 5/5 seeds, RE first-spike SD 0.03 ms.

Then try the levers **one at a time**, keeping each that helps:

1. [x] **Weaker RE gap junctions: `g_gap` 0.03 → ~0.001 µS** (sweep 0.003,
       0.001, 0.0005, 0). Measured cause of RE's microsecond synchrony
       (`Edu-questions.md` Q7, `neuron/re_gap_sync_test.py`): at 0.03 µS one
       junction is 4.3× an RE cell's own conductance, coupling is 0.43 to the
       nearest neighbour and 0.27 across the whole ring, and RE fires within
       0.08 ms even when its input is spread over 25 ms (9.9 ms with gaps
       off). 0.001 µS gives realistic coupling: 0.09 nearest neighbour, 0.008
       across the ring. Check first that RE's spread now follows its input,
       then whether the loop gets a second cycle.
       **Also change it in production:** `g_gap` in `ctx_thalamus_mpi.py`
       (constructor default 0.03). Keep the old value reachable by flag so
       earlier runs can be reproduced.
2. [ ] **Local, topographic RE ↔ TC wiring** instead of all-to-all. Cells on a
       line or ring; each TC gets input from the ~10–20 nearest RE cells, each
       RE from nearby TC cells. Sweep the footprint (5, 10, 20, 40, all).
       Expected effect: phase dispersion, TC firing every 2nd–3rd cycle.
3. [ ] **Cooperative GABA_B on RE → TC**, alongside GABA_A. Add `gabab.mod`
       (Destexhe & Sejnowski 1995 kinetic model; a published version is on
       ModelDB). Check first that one RE spike gives almost no GABA_B current
       and a burst gives a large one. Sweep its conductance from small values
       up; watch that TC isn't silenced (the linear version's failure).
4. [ ] **Ca²⁺-dependent I_h** (`gh_tc`) at physiological levels (around
       Destexhe's `2×10⁻⁵`), for waning and a refractory period of several
       seconds, not for setting the frequency.
5. [ ] Only if still short: RE → TC GABA_A decay (`tau2_re_tc`; an earlier
       sweep showed a modest gain) and RE → RE inhibition strength.

### Results so far (2026-09-24)

Raw output: `res/2026-09-24/ring_gap_sweep.txt`, `ring_loop_sweep.txt`;
figure `out/stage1_ring_thalamus_lfp.png`. `g_gap` is now a CLI flag
(`--g-gap`, default still 0.03); production is **not** changed, since
lever 1 alone did not help.

- **Lever 1 (`g_gap` 0.03, 0.003, 0.001, 0.0005, 0), L6 and tone kicks, 5
  seeds each (50 runs): 1 cycle in every run.** Gap strength does set RE
  synchrony, as Q7 predicted (RE first-spike SD 0.03 → 0.39 ms with an L6
  kick, 0.04 → 0.97 ms with a tone), but RE stays sub-millisecond even at
  g_gap = 0 because every RE cell sums the same common input (100 L6 inputs
  or ~91 of 346 TC cells), and a second cycle never appears.
- **Where the loop breaks (single run, g_gap 0.001):** the kick fires all
  TC, RE bursts (7 spikes / 15 ms), TC is pulled to −83 mV and **does
  rebound ~80 ms later** (spindle timing), but weakly: I_T availability h_T
  only reaches 0.19 from 0.03 at TC's −70 mV rest during the 8 ms GABA_A
  IPSP, so each TC cell fires ~1 spike and the rebound spreads over ~100 ms
  (107, 104, 75, 44 … TC cells per 20 ms). That never recruits RE, which sits
  at −77 mV after its burst.
- **Exploratory (lever 5 early, plus TC → RE gain):** `tau2_re_tc` 8/20/40/80
  ms × `g_tc_re` 0.011/0.03/0.06, g_gap 0.001, 3 seeds: 1 cycle everywhere
  except `tau2_re_tc` = 80 ms with `g_tc_re` = 0.06, which rings for the
  whole 3 s — but at **~4.3 Hz** (RE volleys 200–240 ms apart), the
  GABA_B-like / absence-like thalamic rhythm (von Krosigk et al. 1993), not
  a spindle. Its 10–15 Hz band nevertheless looks like a continuous spindle
  (bottom of the figure): the Stage 0 measure is what exposes it.

**Implication for the next levers.** The missing piece is a **strong,
synchronous TC rebound at ~70–80 ms**: brief (GABA_A-like) inhibition gives
the right timing but too little I_T de-inactivation from a −70 mV rest; long
inhibition gives enough de-inactivation but the wrong (4 Hz) timing. So
cooperative GABA_B (lever 3) is likely to push toward 3–4 Hz, not 12 Hz,
and should be tested with that expectation. Two levers the list above does
not have, both physiological, are worth trying first:
(a) **deeper TC / RE resting hyperpolarisation** (sleep depth — the K⁺ leak
knob from `PLAN-auditory-input.md`, `kl_scale` > 1), so a brief IPSP finds
more I_T available; (b) **RE → TC IPSP amplitude** (`g_re_tc`) at the
8 ms decay, i.e. stronger rather than longer inhibition. Then lever 2
(local wiring) for waxing/waning and TC participation every 2nd–3rd cycle.

**Done when:** a single kick produces, on ≥ 4 of 5 seeds, a train of ≥ 6
cycles at 10–15 Hz lasting ≥ 0.5 s with a growing-then-fading envelope, and a
second kick within ~2 s produces a weaker or no spindle (refractoriness).

## Stage 2 — reconnect the cortex (local, reduced scale)

Only after Stage 1 passes. Use `ctx_thalamus_mpi.py` at `--scale 0.05`–`0.1`.

- [ ] Port the Stage 1 wiring/mechanisms into `ParallelCorticoThalamicNet`
      behind flags (default off, so existing runs are unchanged).
- [ ] Slow the cortical slow oscillation to < 1 Hz (L5 IB pacemaker /
      L5 recurrence parameters).
- [ ] Make the L6 → TC / RE kick weaker or more spread in time (lower
      `g_l6_tc`, `g_l6_re`, or delay jitter), so it triggers a spindle rather
      than resetting the thalamus.
- [ ] Check spindles start in UP states, and that not every UP state carries
      one.

**Done when:** Stage 1 criteria hold with the cortex attached, SO < 1 Hz, and
spindles are nested in UP states.

## Stage 3 — full-scale confirmation on MN5

- [ ] One 200 s run at `--scale 1.65` with the Stage 2 settings.
- [ ] Evaluate with the Stage 0 measure and the existing figures
      (`ctx_analyze.py`, `volley_stats.py`).
- [ ] Update `neuron/README.md` and `Edu-questions.md`.

**Done when:** Stage 1 criteria hold at 5031 cells over 200 s without drift.

## Out of scope for now

- Changing the LFP proxy or computing EEG (`Edu-questions.md` Q3–Q4).
- Matrix (TC → L1) pathway, SST interneurons, dendritic disinhibition.
- Two-threshold spindle detection on the LFP proxy. Stage 0 replaces it for
  judging spindles.

## Risks

- **Stage 1 may still fail.** If local wiring + cooperative GABA_B + I_h don't
  ring, the single-compartment TC/RE cells may be the limit: RE spindle
  bursting is dendritic (Ca_v3.3). The next step then would be a dendritic
  compartment for RE.
- **Earlier results may shift.** The Option 1 / Option 2 choice was scored
  with `n_spindle_events` / `spindle_cv`, which count single volleys
  (`Edu-questions.md` Q5). Re-score with the Stage 0 measure once spindles
  exist.

## References

- Destexhe A, McCormick DA, Sejnowski TJ (1993). A model for 8–10 Hz spindling
  in interconnected thalamic relay and reticularis neurons. *Biophys J*.
- von Krosigk M, Bal T, McCormick DA (1993). Cellular mechanisms of a
  synchronized oscillation in the thalamus. *Science*.
- Destexhe A, Sejnowski TJ (1995). G protein activation kinetics and spillover
  of GABA may account for differences between inhibitory responses in the
  hippocampus and thalamus. *PNAS*.
- Destexhe A, Bal T, McCormick DA, Sejnowski TJ (1996). Ionic mechanisms
  underlying synchronized oscillations and propagating waves in a model of
  ferret thalamic slices. *J Neurophysiol*.
- Bal T, McCormick DA (1996). What stops synchronized thalamocortical
  oscillations? *Neuron*.
- Lüthi A, McCormick DA (1998). Periodicity of thalamic synchronized
  oscillations: the role of Ca²⁺-mediated upregulation of I_h. *Neuron*.
- Bazhenov M, Timofeev I, Steriade M, Sejnowski TJ (2002). Model of
  thalamocortical slow-wave sleep oscillations and transitions to activated
  states. *J Neurosci*.
