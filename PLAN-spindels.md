# Plan: make the thalamus produce real sleep spindles

Status: **Stage 0 done. Stage 1: the isolated thalamus produces
spindle-like trains in the 10–15 Hz band with partial refractoriness
(2026-09-28, strong I_h locking; see "Option 1"). Stage 2 done (2026-09-30):
with the cortex attached, SO 0.93 Hz, spindles nested in UP states, 10–15 Hz,
median 0.55 s, half of them waxing/waning, refractory (see "Stage 2
results" and "Stage 2 follow-up"). Stage 3 first MN5 runs (2026-10-01):
every spindle meets the criteria at 5031 cells over 200 s with no drift, but
the full-scale SO is clockwork and every UP state carries an identical
spindle (see "Stage 3 results"). Round 2 (2026-10-02): with cortical
background input at 200 Hz the SO is irregular (0.77 Hz, CV 0.38), 63% of
UP states carry a spindle, refractoriness is visible; spindles are ~1
cycle short (median 498 ms) (see "Stage 3 round 2").** See "Results so far", "RE recovery" and "Lever 4 (I_h)" in
Stage 1. Background and evidence:
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
2. [x] **Local, topographic RE ↔ TC wiring** instead of all-to-all. Cells on a
       line or ring; each TC gets input from the ~10–20 nearest RE cells, each
       RE from nearby TC cells. Sweep the footprint (5, 10, 20, 40, all).
       Expected effect: phase dispersion, TC firing every 2nd–3rd cycle.
3. [ ] **Cooperative GABA_B on RE → TC**, alongside GABA_A. Add `gabab.mod`
       (Destexhe & Sejnowski 1995 kinetic model; a published version is on
       ModelDB). Check first that one RE spike gives almost no GABA_B current
       and a burst gives a large one. Sweep its conductance from small values
       up; watch that TC isn't silenced (the linear version's failure).
4. [x] **Ca²⁺-dependent I_h** (`gh_tc`) at physiological levels (around
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

### Levers (a) sleep depth and (b) fast-inhibition strength (2026-09-24)

Raw output: `res/2026-09-24/ring_depth_sweep.txt`, `ring_gretc_sweep.txt`,
`ring_combo_sweep.txt`; figure `out/stage1_levers_ab_thalamus_lfp.png`.
All at g_gap 0.001, L6 kick. Sleep depth is the K⁺ leak scale
(`kl_scale_tc` / `kl_scale_re`, new constructor options; 1 = the nrem
preset, TC rest −75.6 → −84.5 mV and I_T availability 0.14 → 0.59 from
1 to 2).

- **(a) Deeper TC / RE rest** (`kl_scale_tc` 1, 1.25, 1.5, 2 ×
  `kl_scale_re` 1, 1.5; 3 seeds): **1 volley in all 24 runs**, and the TC
  rebound **disappears** (10–15 TC cells instead of ~100 per 20 ms at 1.5;
  none at 2). At −81 to −85 mV the RE IPSP (E_GABA −85 mV) barely
  hyperpolarises TC further, and afterwards TC just returns to a deep rest:
  there is nothing to rebound to. Deeper is not better; the current depth
  is near the useful range.
- **(b) Stronger fast inhibition** (`g_re_tc` 0.015, 0.03, 0.06, 0.12 at
  τ 8 ms; 3 seeds): **1 volley in all 12 runs.** The rebound grows only a
  little (peak 104 → 164 TC cells per 20 ms at 8×), since TC is already
  pulled close to E_GABA, and stays spread over ~100 ms; RE is never
  recruited.
- **Combined** (`g_re_tc` 0.015/0.06/0.12 × `g_tc_re` 0.03/0.06/0.12 ×
  τ 8/20 ms; 2 seeds, 36 runs): 1 volley everywhere except the extreme
  corner (all 8×, τ 20 ms), which is **bistable**: seed 0 starts oscillating
  from the start-up transient at 0.12 s and runs at a steady **~8 Hz**
  (RE volleys 119–129 ms apart, TC participation ~32% per cycle) for the
  whole run, untouched by the kick; seed 1 fires once after the kick and
  stays silent. Not a spindle (not triggered, no waxing/waning, never
  terminates), but it is the first regime with TC firing on ~1 cycle in 3,
  and termination is what lever 4 (Ca²⁺-dependent I_h) is for.

**Where this leaves Stage 1.** Five levers now fail the same way: one RE
volley, then a TC rebound that is too weak and too dispersed (~100 ms) to
recruit RE again. RE in this model recovers slowly after its burst (sits at
−77 mV, above its −81 mV rest, for hundreds of ms), and every RE cell sees
the same all-to-all TC input. What is left in the list, in the order the
evidence now suggests:
1. **Lever 2, local RE ↔ TC wiring** — the only lever aimed at the
   dispersion itself (different TC subsets rebounding on different cycles).
2. **RE recovery after its burst** (not in the list): why RE stays
   depolarised; SK2 / I_T2 kinetics or RE → RE inhibition as the RE
   hyperpolarising drive that re-primes RE I_T between cycles.
3. **Lever 4 (I_h)** only for termination, once a regime rings.
4. **Lever 3 (GABA_B)** last: the long-IPSP runs here ring at 4 Hz, the
   known GABA_B / absence rhythm, so it is expected to push the wrong way.
If these fail too, the risk named below applies: the single-compartment RE
cell may be the limit (dendritic Ca_v3.3).

### Lever 2 (local wiring) and RE recovery: the thalamus rings (2026-09-28)

Raw output: `res/2026-09-28/ring_*_sweep.txt`; figures
`out/stage1_local_wiring_lfp.png`, `out/stage1_ringing_thalamus_lfp.png`,
`out/stage1_localkick_kick2_lfp.png`. New options, all defaulting to the
production behaviour: `thal_footprint` (`--thal-footprint`: each TC hears
its F nearest RE on a ring, each RE the TC cells covering the same stretch),
`ek_tc` / `ek_re`, `taur_re`; `thal_ring_test.py` gained `--kick-frac`
(local kick) and `--kick2-t` (second kick).

**Lever 2 alone: negative.** F = 5, 10, 20, 40 × global or local (20% of
ring) kick, 3 seeds: 1 RE volley in all 27 runs. With a local kick the TC
rebound stays in the kicked patch and nothing spreads; with random wiring
the kicked RE patch inhibits TC around the whole ring, which all rebounds,
but RE never fires again.

**The cause: RE cannot burst again for > 200 ms.** Single-cell test (a
second AMPA input Δ ms after a burst): a rested RE cell bursts at 0.02 µS;
after a burst it does not burst at Δ = 80–200 ms even at 0.06 µS. Three
reasons, each measured:
1. **E_K was never set: every K⁺ current in the model uses NEURON's default
   ek = −77 mV.** For RE (rest −81.5 mV) SK2 and hh2 K⁺ therefore
   *depolarise* after a burst: RE sits at −77 mV for ~200 ms and its I_T2
   recovers slowly (h 0.27 → 0.37). Destexhe's thalamic models use −95 to
   −100 mV. (Cortex, rest −70 mV, is less affected; not changed here.)
2. **SK2 stays saturated for ~200 ms** (g = gkbar, 20× the leak) because
   the Ca²⁺ pool clears with τ = 80 ms (`cad` default; Destexhe's cadecay
   uses ~5 ms) and peaks 9× above the SK2 K_d. With ek −95 and τ 20 ms
   RE re-bursts at 100–120 ms.
3. **TC → RE is too weak to re-trigger RE even if TC were synchronous:**
   0.011 µS total, below the ~0.02 µS a rested RE cell needs. Mushtaq 2024
   Table 3 has TC → RE as the *strongest* thalamic synapse.

**With the three RE fixes the kick evokes a spindle-like train** (isolated
MN5 thalamus, g_gap 0.001, 5 seeds). `ek_re = −95`, `taur_re = 5`,
`g_tc_re = 0.06`, local wiring:

| setting | cycles | duration | freq | TC part./cycle | waxing/waning |
|---|---|---|---|---|---|
| global kick, F 5 | 7–10 (5/5 ≥ 6) | 0.58–0.87 s | 10.3–10.6 Hz | 32–49% | no (starts at 100%) |
| global kick, F 10 | 8–11 (5/5) | 0.68–1.01 s | 9.9–10.3 Hz | 34–50% | no |
| global kick, F 20 | 11–13 (5/5) | 0.98–1.16 s | 9.7–10.4 Hz | 33–49% | no |
| local kick (20%), F 10 | 7–8 (5/5) | 0.59–0.72 s | 9.7–10.2 Hz | | **5/5** |
| local kick (20%), F 20 | 8–11 (5/5) | 0.70–0.98 s | 9.7–10.2 Hz | | **5/5** (RE 0.37 → 0.63 → 0.15) |

The network is silent before the kick, trains stop on their own, and with
a local kick the activity spreads along the ring (propagating waves,
Destexhe et al. 1996). Controls: `ek` −95 on TC as well turns ~18% of TC
cells into 15–19 Hz pacemakers with no input, so it is applied to RE only;
`taur_re` 20 ms or `g_tc_re` 0.09 ring for the whole 3 s (30+ cycles) and
sometimes start before the kick. Random wiring with the fixes (and ek −95
on both) gave 4 cycles at ~11 Hz.

**Against the "done when":** ≥ 6 cycles on ≥ 4/5 seeds ✅ (F ≥ 10, also with
a local kick); ≥ 0.5 s ✅; growing-then-fading envelope ✅ with a local
start; 10–15 Hz ⚠️ ~10 Hz, the bottom edge (9.7–10.6); **refractoriness ❌**:
a second kick at 2.5 s (~0.7 s after the first spindle ends) evokes a
full-size second spindle on 5/5 seeds.

**Next:** lever 4, Ca²⁺-dependent I_h (`gh_tc`) at physiological levels, for
refractoriness; the README found I_h also raises the loop frequency, which
may move ~10 Hz into the band. Then Stage 2 with these settings: none of
them is in production yet (all behind options).

### Lever 4: Ca²⁺-dependent I_h (2026-09-28)

Raw output: `res/2026-09-28/ring_ih2_*.txt`; figure `out/stage1_ih_lfp.png`.
Setting: the working point above (`ek_re −95`, `taur_re 5`, `g_tc_re 0.06`,
F 10, g_gap 0.001), kicks at 3.0 s and 4.5 s (a 3 s settle, because I_h
takes ~1 s to reach steady state), 5 seeds. New option `depth_tc` (TC Ca²⁺
pool depth, default unchanged).

- **I_h alone depolarises TC:** at Destexhe's 2e-5 S/cm² TC rests at
  −65.6 mV (I_T availability 0.013 instead of 0.14), and the uncompensated
  network oscillates almost continuously (16–55 cycles at 5e-6). So each
  I_h level is **balanced by the K⁺ leak** (`kl_scale_tc` 1.33 / 1.65 /
  2.27 for 5e-6 / 1e-5 / 2e-5), which restores the −75.6 mV rest exactly.
- **The Ca²⁺ dependence needs a thinner Ca²⁺ pool.** With TC's 10 µm pool a
  spindle's I_T influx peaks at 0.17 µM against ihca's 2 µM
  half-activation (regulating factor p1 ≤ 0.04). With 1 µm (Destexhe-like)
  it reaches 1.5 µM, p1 → 1, 90% of I_h locks open and stays 79% locked
  1.5 s later, holding TC ~2 mV more depolarised (single-cell test).
- **Network results (`depth_tc 1`):**

  | I_h (S/cm²) | cycles | duration | frequency | TC part./cycle | 2nd kick 1.5 s later |
  |---|---|---|---|---|---|
  | 0 | 8–12 | 0.68–1.06 s | 10.0–10.4 Hz | 36% | same (8–12) |
  | **5e-6** | **7–8** | **0.48–0.59 s** | **11.7–12.6 Hz** | 44% | same (7–8) |
  | 1e-5 | 6 | 0.37–0.40 s | 12.6–13.5 Hz | 44% | same (6–7) |
  | 2e-5 | 3 | 0.15–0.16 s | 12.3–13.0 Hz | 38% | 3–4 |

  I_h **shortens** the spindle and **raises its frequency into the band**,
  and removes all pre-kick activity (0 spikes). At 5e-6 the thalamic LFP
  waxes over the first 3–4 cycles and wanes, even with a global kick.
- **No refractoriness at any level.** Likely reason: the L6 kick puts
  0.03 µS directly onto every RE cell, above the ~0.02 µS a rested RE cell
  needs, so RE starts the second spindle whatever TC's state; the
  Ca²⁺-locked I_h acts on TC and only shortens trains.

**Against the "done when"** (I_h 5e-6, leak-compensated, 1 µm pool):
≥ 6 cycles ✅ 5/5; 10–15 Hz ✅ (11.7–12.6); ≥ 0.5 s ⚠️ 4/5 seeds by median
(0.48–0.59 s); envelope ✅; refractoriness ❌.

**Next for refractoriness:** test with a kick that does not force RE — a
**tone kick** (reaches RE only via TC, so TC's I_h state should gate it;
also the case the auditory plan cares about) or an L6 kick near RE's
threshold — and measure the 2nd-spindle probability vs interval (0.5–8 s).
If still absent, add the RE side of refractoriness (review §V.D: RE
Ca²⁺/Na⁺-dependent K⁺ hyperpolarisation).

### Refractoriness with a tone kick (2026-09-28)

Raw output: `res/2026-09-28/ring_tone_refr.*.txt`; figure
`out/stage1_tone_refractoriness_lfp.png`. Same working point, tone kick
(60 dB, 50 ms, IC → TC input; RE reached only via TC) at 3.0 s, a second
identical tone Δ later, 5 seeds; without I_h and with I_h 5e-6
(compensated, 1 µm pool).

- **One tone evokes a spindle** in the isolated thalamus: 7–11 cycles,
  0.68–1.08 s, 8.8–9.6 Hz without I_h; 5 cycles with I_h.
- **No refractoriness.** Second spindle (cycles) vs the first:

  | Δ | no I_h (first: 7–11) | I_h 5e-6 (first: 5) |
  |---|---|---|
  | 0.75 s | 0 — tone 2 arrives during the last cycle of spindle 1, which it ends | 4–5 (tone 0.35 s after the end) |
  | 1.5 s | 6–9 | 5–6 |
  | 3 s | 6–9 | 5–6 |
  | 6 s | 6–9 | 5–6 |

  The TC response to the tone itself is unchanged (3.2 vs 3.2–3.5 spikes/
  cell in 0–50 ms without I_h; 2.5 vs 2.4–2.6 with). The Ca²⁺-locked I_h
  adds only ~2 mV of TC depolarisation, and the tone (like the L6 kick)
  drives TC well past that.

**So refractoriness is the one Stage 1 criterion still missing, and it
does not come from I_h at this strength.** Options, roughly in order:
1. **Stronger / longer I_h up-regulation** at the same resting balance:
   `ginc` (locked-open conductance ratio, 2 now), `k2` (unbinding, τ 2.5 s),
   so the post-spindle depolarisation is several mV and lasts seconds.
2. **RE-side refractoriness** (review §V.D): a slow Ca²⁺- or Na⁺-activated
   K⁺ current in RE that accumulates over a spindle.
3. **Accept it as a network property** and move to Stage 2: in vivo the
   cortex and neuromodulation also gate spindle timing; with the cortex
   attached the 5–10 s spacing may come from the slow oscillation.
None of the Stage 1 settings is in production yet.

### Option 1: stronger I_h locking — partial refractoriness (2026-09-28)

Raw output: `res/2026-09-28/ring_ihlock.*.txt`; figure
`out/stage1_ih_locking_lfp.png`. New options `taur_tc` (TC Ca²⁺ clearance)
and `ginc_tc` (locked-open conductance ratio); defaults unchanged.

- **Why the previous I_h runs could not be refractory:** with the 1 µm pool
  and TC's 80 ms Ca²⁺ clearance, window I_T holds resting Ca²⁺ at 1.1 µM, so
  I_h is ~83% locked *at rest* and a spindle cannot lock it further
  (`k2` had no effect for the same reason). With Destexhe's ~5 ms clearance,
  26% locked at rest → 77% after a spindle.
- **Single-cell screen** (τ 5 ms, rest held at −75.6 mV by the K⁺ leak):
  after a 7-cycle 12 Hz RE barrage TC stays depolarised by +2.9 / +5.2 /
  +7.8 mV for > 4 s at `ginc` 2 / 4 / 8 (`kl_scale_tc` 1.273 / 1.422 /
  1.671).
- **Network, tone kick, 3 seeds** (`gh_tc 5e-6`, `depth_tc 1`, `taur_tc 5`,
  plus the working point):

  | ginc | 1st spindle | 2nd tone +1.5 s | 2nd tone +4 s |
  |---|---|---|---|
  | 2 | 5–7 cyc, 0.41–0.60 s, ~10 Hz | 6–7 cyc (no change) | 6–7 cyc |
  | 4 | 7–9 cyc, 0.51–0.67 s, ~11.7 Hz | 7–9 cyc, ~10% shorter | 8–9 cyc |
  | **8** | **8 cyc, 0.49–0.55 s, 12.7–14.3 Hz** | **6 cyc, 0.45 s, ~11 Hz** | **5–6 cyc, 0.43–0.44 s** |

  At `ginc` 8 a second tone 1.5–4 s later gives a spindle with 25–35% fewer
  cycles in 3/3 seeds: a **weaker** second spindle, which meets the Stage 1
  wording ("weaker or no spindle"), but it is never abolished. Stronger
  locking also speeds the first spindle up (volley interval 94 → 78 →
  65 ms at `ginc` 2 → 4 → 8).

**Stage 1 status with `ginc` 8:** ≥ 6 cycles ✅ (8, 3/3); 10–15 Hz ✅
(12.7–14.3); ≥ 0.5 s ⚠️ (0.49–0.55); waxing/waning ✅ (earlier, local
start); refractoriness ✅ partial (weaker, not absent). Only 3 seeds and
one level (60 dB); the full-refractoriness question (no spindle within
~2–5 s) is still open — RE-side K⁺ accumulation (option 2) is the
remaining lever. Next logical step: Stage 2 (cortex attached) with these
settings, all still behind options.

**Done when:** a single kick produces, on ≥ 4 of 5 seeds, a train of ≥ 6
cycles at 10–15 Hz lasting ≥ 0.5 s with a growing-then-fading envelope, and a
second kick within ~2 s produces a weaker or no spindle (refractoriness).

## Stage 2 — reconnect the cortex (local, reduced scale)

Only after Stage 1 passes. Use `ctx_thalamus_mpi.py` at `--scale 0.05`–`0.1`.

- [x] Port the Stage 1 wiring/mechanisms into `ParallelCorticoThalamicNet`
      behind flags (default off, so existing runs are unchanged).
- [x] Slow the cortical slow oscillation to < 1 Hz (L5 IB pacemaker /
      L5 recurrence parameters).
- [x] Make the L6 → TC / RE kick weaker or more spread in time (lower
      `g_l6_tc`, `g_l6_re`, or delay jitter), so it triggers a spindle rather
      than resetting the thalamus.
- [x] Check spindles start in UP states, and that not every UP state carries
      one.

**Done when:** Stage 1 criteria hold with the cortex attached, SO < 1 Hz, and
spindles are nested in UP states.

### Stage 2 results (2026-09-29)

Raw output: `res/2026-09-29/stage2/stage2_report_*.txt`; figures
`out/stage2_so_spindles_lfp.png`, `out/stage2_l6_kick_lfp.png`. Measure:
`neuron/stage2_report.py` (UP-state onsets from the L5E + L6E rate; Stage 0
RE volley trains; spindle = ≥ 6 cycles; "in an UP state" = starts −50 …
+300 ms from an onset; `--pool` over seeds).

Setup: `ctx_thalamus_mpi.py --scale 0.1 --thal-scale 1.65` (MN5-size
thalamus, 346 TC / 91 RE, under a 0.1-scale Option 2 cortex: `--g-i-e-l5
0.02 --g-l5-rec 0.013 --tau2-l5-rec 200 --taur-l5e-rs 120 --l5-rec-mech
nmda`), 20 s runs, first 2 s skipped. The Stage 1 thalamus is `--g-gap 0.001
--thal-footprint 10 --ek-re -95 --taur-re 5 --g-tc-re 0.06 --gh-tc 5e-6
--kl-scale-tc 1.671 --depth-tc 1 --taur-tc 5 --ginc-tc 8`. New options (all
default off): `--thal-scale`, `--taur-l5-ib`, `--het-seed`,
`--l6-delay-spread`, and the Stage 1 flags.

1. **The Stage 1 thalamus survives the cortex.** Production thalamus: 62 RE
   volleys, 61 of them single (0 spindles). Stage 1 thalamus: 8 spindles
   at 10–12 Hz, but up to 2.5 s long (the production L6 kick 0.03 keeps
   re-driving it). A weaker L6 kick (0.01 to TC and RE) brings them down
   to 0.4–0.9 s.
2. **SO < 1 Hz via the L5 IB Ca²⁺ clearance** (`--taur-l5-ib`, default
   500 ms):

   | taur_l5_ib | SO (Hz) |
   |---|---|
   | 500 | 1.89 |
   | 1000 | 1.28 |
   | 1500 | 1.11 |
   | **2000** | **0.91 pooled over 5 seeds (0.78–1.11)** |
   | 3000 | 0.56 (few spindles) |

   UP-state intervals are irregular (CV 0.6–0.9).
3. **Nesting ✅ and refractoriness ✅** (2000 ms, L6 0.01, 5 seeds, 90 s):
   - 100% of spindles start in an UP state;
   - 39% of 82 UP states carry one (not every UP state);
   - an UP state right after a spindle-carrying one carries a spindle 27% of
     the time, vs 47% after one without; the median spindle-to-spindle
     interval is 2.5 s.
4. **The envelope problem.** The 0.1-scale cortex's UP onset is a
   population spike: all 42 L6E cells fire within ~2 ms. With L6 → TC and
   RE both at 0.01 it fires every TC and RE cell in the first cycle (4.2 RE
   spikes/cell, 0 ms spread), and each following cycle is smaller: the
   spindles are **decrementing**, not waxing (waxing/waning in 9% of 32).
   Variants at 2000 ms:

   | L6 → TC / RE | seeds | spindles | ≥ 0.5 s | median | wax/wane | TC part./cycle |
   |---|---|---|---|---|---|---|
   | 0.01 / 0.01 | 5 | 32 | **81%** | 583 ms | 9% | 68% |
   | 0.01 / 0.01, + 0–40 ms delay spread | 1 | 5 | 60% | 537 ms | 0% | 59% |
   | 0.005 / 0.005 | 1 | 4 | 25% | 480 ms | 0% | 71% |
   | 0.005 / 0.005, + spread 40 | 1 | 1 (15.9 Hz) | 100% | 503 ms | 0% | 59% |
   | 0.003 / 0.003, + spread 40 | 1 | 0 (no RE volleys) | – | – | – | – |
   | 0.005 / 0.01 | 1 | 3 | 33% | 442 ms | 0% | 65% |
   | **0.003 / 0.01** | **4** | **10** | 10% | 400 ms | **50%** | **49%** |
   | 0.003 / 0.015 | 1 | 1 | 0% | 339 ms | 0% | 25% |

   - Spreading the L6 arrival in time does not help.
   - Making the cortex drive **mainly RE** (L6 → TC 0.003, L6 → RE 0.01)
     gives the textbook sequence (review: cortex → RE → TC). RE fires first
     and TC is almost silent in cycle 1 (0–25%). TC then rebounds with
     ~50% participation per cycle, i.e. each TC cell fires about every 2nd
     cycle, as the criteria ask. The RE and LFP envelope grows over 2–3
     cycles and then fades.
   - The cost is **shorter spindles**: 5–7 cycles, ~0.4 s, and only 18% of
     UP states carry a ≥ 6-cycle one.

**Against the "done when":**
- SO < 1 Hz ✅
- nested in UP states ✅ (100%), not every UP state ✅
- refractoriness ✅ (spontaneous, partial)
- 10–15 Hz ✅ (84–90%)
- Stage 1 criteria ⚠️: either ≥ 0.5 s (81%, L6 0.01 / 0.01) or
  waxing/waning plus TC every 2nd cycle (50% / 49%, L6 0.003 / 0.01), not
  both in one setting.

The spindle length is capped by the isolated thalamus's own ~0.5 s (Stage 1
`ginc` 8). The RE-first start spends one of those cycles, which pushes it
under 0.5 s.

**Next (small, before Stage 3):** lengthen the intrinsic spindle while
keeping the RE-weighted L6 input:
- `ginc_tc` 4 (Stage 1: 7–9 cycles, 0.51–0.67 s, still partly refractory);
- or I_h 5e-6 → 3e-6 (fewer, longer spindles);
- then 5 seeds.

A larger cortex may also help, if its UP onsets are less synchronous; that
is testable in Stage 3 at `--scale 1.65`.

### Stage 2 follow-up: background input and weaker I_h locking (2026-09-30)

Raw output: `res/2026-09-30/stage2/stage2_report_*.txt`; figure
`out/stage2_ginc4_noise_lfp.png`. Same setup as above with `--taur-l5-ib
2000`. New options (default off): `--cx-noise-rate` / `--cx-noise-w` and
`--het-seed`. `stage2_report.py` now also prints the UP-onset spread (SD and
10–90% range of L6E first spikes).

**How the SO is generated** (checked on single cells):
- **No external drive.** There is no periodic input, noise or current
  injection.
- **The pacemaker is the L5 IB cell.** Its persistent Na⁺ current drives a
  burst (the UP state), SK2 ends it, and the slow Ca²⁺ clearance
  (`taur_l5_ib`) sets the DOWN state. A lone IB cell with no input bursts
  every 1.24 / 2.4 / 4.7 s at 500 / 1000 / 2000 ms.
- **The network speeds this up 2–4×.** The first IB cell to recover recruits
  the rest through the IB gap junctions and the L5 NMDA recurrence.
- **The thalamus follows and does not set the pace.** The SO is the same
  with the production and Stage 1 thalamus.
- **Consequence:** the DOWN → UP transition is deterministic, and the column
  rises within ~2 ms (L6E 10–90% range 1.5–2.9 ms).

1. **Background input does not desynchronise the UP onset.** A Poisson AMPA
   train goes to every cortical E cell: 2e-4 µS, a 0.7 mV EPSP; at
   100–200 Hz it gives 1–1.3 mV of membrane noise and a lone cell stays
   silent.
   - The onset range grows only from ~2 to ~3.5 ms. The gap junctions and
     recurrence still recruit the column at once.
   - It ends DOWN states earlier: SO 1.1–1.3 Hz at IB 2000 ms, and
     **1.6 Hz at IB 3000 ms**. Once the input is present, the SO period is
     set by noise-triggered recurrence, not by the IB clearance.
   - Some spindles now start in part of the RE ring and spread (visible in
     the figure).
   - Not adopted: it does not fix the envelope, and it breaks SO < 1 Hz.
2. **Weaker I_h locking restores the length.** `--ginc-tc 4
   --kl-scale-tc 1.422` plus the RE-weighted cortical input (`--g-l6-tc
   0.003 --g-l6-re 0.01`), 5 seeds, 90 s:

   | Criterion | Result |
   |---|---|
   | SO | **0.93 Hz** pooled (0.72–1.17 per seed; 1 of 5 seeds > 1 Hz) |
   | Spindles start in an UP state | 88–100% per seed |
   | UP states carrying a spindle | 43% of 84 |
   | ≥ 6 cycles | 37 spindles, 6–10 cycles |
   | 10–15 Hz | 89% (a few at 15.5–15.8 Hz) |
   | ≥ 0.5 s | 73%, median 545 ms; seed medians ≥ 0.5 s in 4/5 seeds |
   | Waxing/waning | 51% (seeds: 83 / 0 / 78 / 25 / 62%) |
   | TC participation per cycle | 51% (each TC cell fires about every 2nd cycle) |
   | Refractoriness | P(spindle \| previous UP had one) 19% vs 58% otherwise; min interval 0.85 s, median 2.7 s |

   Compared with `ginc` 8, the spindle loses less of its drive to I_h, so it
   gains 1–2 cycles. That makes up for the cycle the RE-first start costs.

**Stage 2 status: done.** Every "done when" criterion holds pooled over 5
seeds:
- SO < 1 Hz;
- spindles nested in UP states, not every UP state carrying one;
- Stage 1 criteria, with waxing/waning in half of the spindles and not all
  seeds, and the length criterion met in 4/5 seeds.

Settings for Stage 3:

    --taur-l5-ib 2000 --g-gap 0.001 --thal-footprint 10 --ek-re -95
    --taur-re 5 --g-tc-re 0.06 --gh-tc 5e-6 --kl-scale-tc 1.422 --depth-tc 1
    --taur-tc 5 --ginc-tc 4 --g-l6-tc 0.003 --g-l6-re 0.01

plus the Option 2 cortex flags.

**Caveats for Stage 3:**
- At `--scale 1.65` the L6 → thalamus convergence and the IB population size
  change, so the L6 weights and the SO rate need checking.
- The 0.1-scale cortex has only 42 L6E cells.
- The synchronous UP onset is a property of the cortical wiring (IB gap
  junctions, all-to-all-like recurrence), not of the thalamus.
  Topographic cortical wiring is the lever if propagating UP states are
  wanted.

## Stage 3 — full-scale confirmation on MN5

- [ ] One 200 s run at `--scale 1.65` with the Stage 2 settings.
      Prepared 2026-09-30: `run_ctx_nrn.sh SPINDLE=stage2` (commands in
      `MN5_NEURON.md`, "Spindle Stage 3"). A 12 s local preview at full
      scale keeps the spindles: 8 cycles, 13.6 Hz, 0.51 s, 75%
      waxing/waning. But the SO slows to 0.40 Hz and turns clockwork (every
      UP state carries a spindle), so the MN5 submission adds
      `--taur-l5-ib` 1500 and 1000 variants.
- [ ] Evaluate with the Stage 0 measure and the existing figures
      (`ctx_analyze.py`, `volley_stats.py`).
- [ ] Update `neuron/README.md` and `Edu-questions.md`.

**Done when:** Stage 1 criteria hold at 5031 cells over 200 s without drift.

### Stage 3 results: first MN5 runs (2026-10-01)

Raw output: `res/2026-10-01/s3/stage3_report.txt`, `stage3_drift.txt` (the
`.npz` files are local only); figure `out/stage3_mn5_lfp.png`. Setup:
- 5031 cells (`SCALE=1.65`), 200 s, 100 ranks on MN5, ~10 min wall each;
- Option 2 cortex + `SPINDLE=stage2`;
- three L5 IB clearances (`taur_l5_ib`).

| | IB 2000 (preset) | IB 1500 | IB 1000 |
|---|---|---|---|
| SO | 0.35 Hz | 0.45 Hz | **0.70 Hz** |
| UP-state interval CV | 0.11 | 0.01 | 0.01 |
| Spindles (≥ 6 cycles) | 68 (20.6/min) | 90 (27.3/min) | 138 (41.8/min) |
| Cycles | 7–9 (median 8) | 7–9 (8) | 7–9 (8–9) |
| Frequency | 13.9 Hz (99% in 10–15 Hz) | 13.9 Hz (98%) | 14.0 Hz (100%) |
| ≥ 0.5 s | 97% (median 506 ms) | 96% (504 ms) | 96% (528 ms) |
| Waxing/waning | 99% | 99% | 99% |
| TC participation per cycle | 50% | 49% | 49% |
| Start in an UP state | 100% | 100% | 100% |
| UP states carrying a spindle | **99%** | **100%** | **99%** |

**Drift:** none. The 40 s windows keep the same UP rate, spindle count,
cycles, frequency and duration from 0 to 200 s, and population rates stay
within ±5% (TC 1.5–3.2, RE 5.8–12.1, L5E 1.3–2.8 Hz/cell).

**What holds:** every individual spindle meets the Stage 1 criteria at full
scale, more cleanly than at scale 0.1:
- 7–9 cycles at ~14 Hz, ~0.5 s;
- waxing/waning in 99%;
- each TC cell fires about every 2nd cycle;
- always nested in an UP state, starting at the DOWN → UP transition.

**What fails:**
1. **The SO is clockwork:** UP-state interval CV 0.01–0.11, against
   ~0.3–0.6 in real slow-wave sleep. At scale 0.1 it was 0.6–0.9; that
   irregularity came from the small network, not from a mechanism.
2. **Every UP state carries a spindle, and they are all the same.** So
   "not every UP state carries one" fails, and **refractoriness cannot be
   tested**: the shortest interval (1.4 s at IB 1000) is longer than the
   thalamus's recovery, so it is never challenged.
3. **The spindle rate is too high.** 21–42/min, against ~2–7/min in human
   N2/N3. This follows directly from 2.
4. ~14 Hz is at the top of the band (fast spindles, 12–15 Hz, are the
   type coupled to SO UP states, so this is acceptable).

**Cause:** the full-scale cortex is a deterministic, fully synchronised
oscillator:
- 437 L5 IB cells, coupled by gap junctions;
- L6E UP-onset spread 1.2 ms (10–90%) in all three runs.

So the thalamus gets an identical kick every cycle, from an identical
recovered state, and answers with an identical spindle. The thalamus is no
longer the problem; the variability has to come from the cortex.

**Next (MN5, cheap, ~10 min per run):** make the full-scale SO irregular.
- **Background input to cortical E cells** (`--cx-noise-rate 100 / 200`).
  At scale 0.1 it sped the SO up and made it irregular. At full scale the
  IB 2000 SO is 0.35 Hz, so speeding it up is welcome.
- **Weaker IB gap junctions** (`--g-l5-gap`, new flag, default 0.02), so
  the IB population can drift partly out of phase.
- **Target:** SO 0.5–1 Hz with CV ≥ 0.3, 20–60% of UP states carrying a
  spindle, and refractoriness measurable as at scale 0.1 (Stage 2:
  19% vs 58%).

### Stage 3 round 2: irregular SO at full scale (2026-10-02)

Raw output: `res/2026-10-02/s3b/stage3b_report.txt`,
`stage3b_drift_n200.txt` and the SLURM logs; figure
`out/stage3b_noise200_lfp.png`. Same setup as round 1 (5031 cells, 200 s,
100 ranks, ~10 min each, `SPINDLE=stage2`, IB 2000 unless noted).

| | preset (ref.) | noise 100 | **noise 200** | IB 1500 + noise 100 | gap 0.005 | gap 0.005 + noise 100 |
|---|---|---|---|---|---|---|
| SO | 0.35 Hz | 0.38 | **0.77** | 0.49 | 0.35 | 0.38 |
| UP interval CV | 0.11 | 0.10 | **0.38** | 0.01 | 0.11 | 0.00 |
| UP states with a spindle | 99% | 99% | **63%** | 100% | 99% | 100% |
| Spindles / min | 20.6 | 22.4 | 29.4 | 30.0 | 20.9 | 22.7 |
| ≥ 0.5 s (median) | 97% (506 ms) | 99% (566) | 45% (498) | 97% (565) | 97% (505) | 100% (566) |
| 10–15 Hz | 99% | 97% | 93% | 96% | 99% | 93% |
| Waxing/waning | 99% | 99% | 98% | 99% | 99% | 99% |
| TC part./cycle | 50% | 47% | 45% | 47% | 50% | 48% |
| P(spindle \| prev. UP had one) vs none | 99% vs – | 99% vs – | **56% vs 75%** | 100% vs – | 99% vs – | 100% vs – |

1. **Weaker IB gap junctions do nothing** (0.005 vs 0.02: identical SO and
   spindles). The full-scale synchrony comes from the recurrent L5 NMDA
   excitation, not the gap junctions.
2. **Background input has a threshold.** At 100 Hz nothing changes. At
   200 Hz the cortex switches regime:
   - DOWN states end stochastically: SO 0.77 Hz with CV 0.38, in the
     range of real slow-wave sleep (~0.3–0.6);
   - spindles **start locally** in part of the TC/RE ring and differ from
     one UP state to the next (see the raster);
   - 37% of UP states carry no spindle;
   - **refractoriness becomes visible**: an UP state right after a
     spindle carries one 56% of the time, vs 75% otherwise.
3. **Drift with noise 200 (40 s windows):** UP rate 0.70–0.82/s,
   18–22 spindles per window, 59–68% of UP states carrying one, rates
   within ±13%, no trend. Duration rises slightly, 460–470 ms in the first
   80 s to ~500 ms after. The SLURM check gives RE events 0.566 vs
   0.614 Hz in the first vs second half (+8%).

**Stage 3 status with `SPINDLE=stage2 --cx-noise-rate 200`:**
- SO < 1 Hz and irregular ✅
- nested in UP states ✅ (100%), not every UP state ✅ (63%)
- refractoriness ✅
- 10–15 Hz ✅ (93%), waxing/waning ✅ (98%), TC every ~2nd cycle ✅ (45%)
- no drift ✅ (fluctuations, no trend)
- **≥ 0.5 s ⚠️:** median 498 ms, 45% ≥ 0.5 s. The background input
  shortens spindles by about one cycle (7–8 instead of 8–9), because they
  start locally.
- **Spindle rate is still high:** 29/min, against ~2–7/min in human
  N2/N3 scalp EEG. That comparison is loose: the model is a single column
  and the scalp counts spindles visible over large areas.

**Next:**
- Repeat noise 200 with 2 seeds (`HET_SEED` 1, 2) to confirm.
- Bracket the noise level (150, 250) for the fraction of UP states
  carrying a spindle.
- Add one cycle of length, e.g. a slightly weaker I_h locking (`ginc`
  3, which needs a local re-balance of `kl_scale_tc` first) or
  `--g-l6-re 0.012`.

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
