# NEURON migration (work in progress)

A staged migration of the auditory thalamo-cortical spindle model from NEST to
NEURON, motivated by [`docs/spindle_review_mapping.md`](../docs/spindle_review_mapping.md)
§5.8: sleep-spindle bursting is a **dendritic** phenomenon (TRN via Ca_v3.3 in
distal dendrites, TC via Ca_v3.1 in primary dendrites; Fernandez & Lüthi 2020
§V), which NEST's point neurons (`ht_neuron`, `aeif_cond_exp`) cannot represent.
That is what stalled in-network burst synchrony. NEURON represents multi-
compartment cells, the T-type Ca²⁺ current, gap junctions, and the extrasynaptic
GABA_A / GABA_B receptor cascade natively.

## Why NEURON (and why not just NESTML)

Both were checked as installable here (NEURON 9.0.1 has a native `cp314` arm64
wheel; NESTML 8.3.0 is pure Python + a C++ build step). The deciding factor is
that **NESTML still generates point neurons** — it would give better channel
kinetics and gap junctions but no dendrites, so it risks landing at the same
wall. NEURON also has a published reference implementation of the Mushtaq et al.
2024 model whose parameters we already mapped. Note NEURON is **not faster** than
NEST for point networks — it buys biophysical fidelity, not throughput; at our
scale both are fast enough.

## Layout

| path | purpose |
|------|---------|
| `mod/itd.mod` | **Destexhe et al. 1996 low-threshold T-current** (GHK + φ-scaled kinetics) — the relay burst mechanism |
| `mod/hh2.mod` | Traub-Miles fast Na⁺ / K⁺; fires repetitively on a plateau (no depol. block) |
| `mod/cav3.mod` | earlier Huguenard-McCormick T-current (ohmic, unscaled τ) — superseded by `itd` |
| `mod/cad.mod` | submembrane Ca²⁺ pool in a private `sk` ion (feeds SK2 **and** `ihca`) |
| `mod/sk2.mod` | SK2 Ca²⁺-activated K⁺ burst terminator |
| `mod/ihca.mod` | **Ca²⁺-dependent I_h** (Destexhe et al. 1996) — the spindle terminator / refractory mechanism |
| `mod/gap.mod` | electrical (gap-junction) coupling between TRN cells (connexin-36) |
| `mod/ical.mod` | HVA Ca²⁺ current (Reuveni et al. 1993 kinetics) — cortical pyramidal spike-triggered Ca²⁺ influx, the SK2 Ca²⁺ source in cortex |
| `mod/inap.mod` | Persistent Na⁺ current — L5 intrinsic bursting (`PYCellIB`), the column's own slow-oscillation pacemaker |
| `tc_neuron.py` | `TCCell` (relay), `RECell` (reticular), `gap_junction()` + demos |
| `tc_network_nrn.py` | `ThalamicNet` — the TC↔RE loop network; emits the `tc_validate` contract |
| `cortex_neuron.py` | `PYCell` (pyramidal, hh2 + ical + cad + sk2 adaptation), `FSCell` (PV⁺ fast-spiking), `CorticalColumn` (L4/L2·3/L5/L6) |
| `ctx_thalamus_network.py` | `CorticoThalamicNet` — cortical column + thalamic TC↔RE loop, closed via TC→L4 and L6→TC/RE |
| `tc_neuron_figures.py` | six-panel spindle figure (`out/tc_neuron_spindles.png`): cell bursts → gap-junction sync → network raster / LFP / spectrogram |
| `arm64/` (git-ignored) | compiled mechanisms, from `nrnivmodl mod` |

## Build & run

```bash
# from the repo root, using the migration venv (see repo .gitignore)
python3 -m venv --system-site-packages .venv-neuron
.venv-neuron/bin/pip install neuron
cd neuron
../.venv-neuron/bin/nrnivmodl mod        # compile the T-current
../.venv-neuron/bin/python tc_neuron.py  # single-cell demo
```

## Status

**Done**
- NEURON toolchain verified end-to-end (a single HH soma fires; NMODL compiles).
- The validator [`tc_validate.py`](../tc_validate.py) is **simulator-agnostic**:
  it scores any backend that emits the documented result contract
  (spikes/traces/meta), imports and runs with **no NEST**, and gives
  byte-identical results under both interpreters (`--self-test`).
- `cav3.mod` implements a correct T-current: it **de-inactivates under
  hyperpolarisation** (h → 1.0) and drives a strong **inward Ca current that
  scales with gCaT** (−0.6 to −1.4 mA/cm²) on release — the exact mechanism the
  NEST point models could not express.
- **Physiological rebound BURST achieved.** A single-compartment cell with
  Traub–Miles `hh2` spikes + the **Destexhe et al. 1996 T-current** (`itd`, GHK
  driving force + temperature-scaled kinetics) fires a genuine thalamic
  low-threshold Ca²⁺ spike:

  | pcabar | I_Ca peak | spikes | freq | burst |
  |---|---|---|---|---|
  | 0 | 0.000 | 0 | — | — (control) |
  | 1.0e-4 | −0.054 | 8 | 257 Hz | 27 ms |
  | **1.7e-4** | **−0.089** | **6** | **161 Hz** | **31 ms** (default) |
  | 2.5e-4 | −0.126 | 3 | 49 Hz | 41 ms |

  The default is a **6-spike burst at 161 Hz over 31 ms** — in the review's 2–6
  spike range, at physiological intra-burst frequency and LTS duration.
  Conductance-based, no depolarisation block, and **not** a phenomenological fit
  (unlike AdEx). Two things were decisive: (1) a **single compartment** (the
  2-compartment version re-primed the T-current via somatic back-propagation,
  giving a ~250 ms plateau); (2) the **φ (temperature) scaling** of the τ's,
  which the earlier `cav3` omitted — without it the LTS was ~3× too long.

- **Reticular (RE/TRN) cell + gap junctions done.** `RECell` (I_T + SK2 + hh2,
  no I_h) fires a **5-spike burst at ~85 Hz** (review: TRN bursts 2 to >10
  spikes). `gap_junction()` electrically couples two RE cells; driving cell A
  alone, the coupling **recruits and synchronises** cell B:

  | g_gap (µS) | A spikes | B spikes | A–B coincidence |
  |---|---|---|---|
  | 0 | 18 | 0 | — (B silent) |
  | 0.002 | 18 | 0 | — |
  | 0.008 | 19 | 7 | 32% |
  | 0.02 | 14 | 14 | **100%** |

  Gap junctions require waveform relaxation, which `ht_neuron` lacks — so this
  TRN synchrony mechanism (Fernandez & Lüthi §V.C.1) is **only** available in the
  NEURON port.

- **Network: 4/5 criteria, both cell types bursting — but the loop runs at 5 Hz, not spindle frequency.** The TC↔RE loop
  (`tc_network_nrn.py`, `ThalamicNet`) now runs the real mechanism in-network —
  up from 2/5:

  | criterion | measured | |
  |---|---|---|
  | RE spikes/burst | **4.35** | ✅ |
  | TC spikes/burst | **3.35** | ✅ |
  | RE V_m < −55 mV | 92% | ✅ |
  | TC V_m < −65 mV | 91% | ✅ |
  | spindle density (events ≥0.5 s) | 0.0/min | ❌ |

  Two fixes got there, both diagnosed by measurement:

  1. **A corticothalamic volley onto TC, not just RE.** The relay had *no
     excitatory input at all*, so once hyperpolarised for I_T de-inactivation it
     could never reach threshold — at `e_pas = −88` it fell silent entirely.
     L6 excites both (≈2:1, Mushtaq Table 3).
  2. **Sleep-state hyperpolarised rest** (`tc_e_pas = −80`, `re_e_pas = −82`).
     The I_T inactivation gate measured **h = 0.05** in-network — 95%
     inactivated, rebound impossible. TC must also stay **above the GABA_A
     reversal (−85 mV)**, or reticular input depolarises instead of inhibiting:
     a narrow window.

- **Kinetics verified against the published model.** `mod/itd.mod` checked
  term-by-term against Destexhe et al. 1996 (`TC.tem` / `IT.mod`):

  | quantity | max relative difference, −100…−40 mV |
  |---|---|
  | m∞ | **0.00%** |
  | h∞ | **0.00%** |
  | τ_m | 2.05% → corrected, now exact |
  | τ_h | 13.9% (>−80 mV branch) → corrected, now exact |

  The compiled mechanism measures τ_h ≈ 87 ms at −85 mV against the formula's
  85.4 ms.

### ⚠️ Correction: the loop oscillates at 5 Hz, not 12 Hz

An earlier version of this file (and PR #37) reported the thalamic loop running
at **12.0 Hz, "in the spindle band"**. That was **wrong** — it came from calling
`detect_peak(rate, fs, 6.0, 20.0)`, i.e. searching a window that *excluded the
true peak*. Measured three independent ways:

| method | result |
|---|---|
| direct population inter-burst interval | 195 ms → **5.1 Hz** (n=91, IQR 173–212 ms) |
| spectrum, unconstrained 1–30 Hz | **5.0 Hz** |
| TC alone / RE alone | 5.0 / 5.1 Hz |
| mean-field V_m spectrum | **5.0 Hz** |

and the reported frequency simply tracks the search window:

| search range | "peak" found |
|---|---|
| 1–30 Hz (honest) | **5.0 Hz** |
| 6–20 Hz | 7.0 Hz |
| 9–16 Hz | 11.0 Hz |

This is the same error class caught earlier in the NEST model (a "9.6 Hz"
detection sitting at the edge of its window) — **always search unconstrained
first**.

**What this means.** The cells burst correctly (that part stands: TC 3.35, RE
4.35 spikes/burst, driven by a genuine Ca_v3.1 current). But the *network*
rhythm is ~5 Hz — the **delta/theta** range, not the 10–15 Hz spindle band. The
period (~195 ms) is close to single-cell I_T de-inactivation (τ_h ≈ 85 ms) plus
burst duration and synaptic delays, which suggests the network is expressing the
**intrinsic I_T/I_h delta-like oscillation** rather than the faster RE↔TC *loop*
rhythm that produces spindles. Destexhe's thalamic models generate both regimes;
distinguishing them is the next real question.

This also fully explains the "duration" failure: at 5 Hz a 0.15 s event is under
one cycle, and the 10–15 Hz detector was reading harmonics of a 5 Hz rhythm.

**Superseded framing — event duration.** Events last ~0.15 s (≈2 cycles
at 12 Hz) against the review's 0.5–3 s, so no event passes the ≥0.5 s criterion
and density reads 0.

**Correction:** an earlier version of this file claimed I_T recovery takes
~320 ms against an 83 ms cycle, making sustained bursting impossible. **That was
wrong** — it omitted the temperature factor `phi_h = 3.0^1.2 = 3.74`:

  | V | τ_h at 36 °C | recovery per 83 ms cycle |
  |---|---|---|
  | −80 mV | 73 ms | 68% toward h∞ = 0.32 |
  | −85 mV | **85 ms** | 62% toward h∞ = 0.62 |
  | −90 mV | 79 ms | 65% toward h∞ = 0.85 |

τ_h is *comparable to* the cycle, so I_T recovers substantially each cycle and
the duration limit has some other, still-unidentified cause.

**Eight hypotheses tested and ruled out** (all leave duration at ~0.15 s):

| # | hypothesis | result |
|---|---|---|
| 1 | I_T recovery kinetics too slow | τ_h ≈ 85 ms ≈ one cycle — not limiting |
| 2 | loop gain (`g_tc_re`, `g_re_tc`) | no effect over 7× range |
| 3 | population size (10→40) | no effect |
| 4 | gap-junction strength (0.03→0) | no effect |
| 5 | heterogeneous resting potentials | no effect |
| 6 | progressive recruitment (0→800 ms) | no effect |
| 7 | ~~**SK2 on RE** (burst terminator, 0→0.03)~~ | ~~no effect~~ — **INVALID, see below** |
| 8 | **GABA_B on RE→TC** (slow inhibition) | max event 0.17→0.25 s, but TC bursting collapses 3.24→0.00 |

GABA_B is implemented and available (`g_re_tc_b`, **default 0**): the review
(§V.B.3) identifies it as the slow inhibition pacing spindle waxing/waning, and
it does lengthen the longest event slightly — but at these conductances it
hyperpolarises TC past what the corticothalamic drive can overcome, so relay
bursting collapses. A narrow useful window may exist between 0 and 0.003.

This is now a genuine open research question rather than an untuned parameter.

## The SK + Ca²⁺-dependent mechanism (Aug 2026)

### ⚠️ Correction: hypothesis #7 was invalid — SK2's Ca²⁺ sensor was never in range

`sk2.mod` had **`kd = 0.5`**, a **µM value left in a mM field** — 1000× too high.
The submembrane pool peaks near 3.7×10⁻² mM, so with `hill = 4`:

| | SK2 open fraction |
|---|---|
| peak activation, old `kd = 0.5` | **3×10⁻⁵** (channel effectively absent) |

**SK2 never opened.** The earlier finding that sweeping `gkbar_sk2` from 0 → 0.03
"had no effect on duration *or* bursting" was therefore an artifact of a dead
channel, not evidence against SK2. Hypothesis #7 is withdrawn.

Two fixes, each verified by measurement:

1. **`kd` → 5×10⁻⁴ mM** (0.5 µM; Hirschberg et al. 1998, Xia et al. 1998).
2. **`depth_cad` → 10** — with `depth = 1` the *resting* pool sat at 6×10⁻⁴ mM
   (600 nM) from the standing T-window current, already above any physiological
   K_d, which pinned SK2 open and abolished bursting. At `depth = 10` resting
   [Ca²⁺] is ~2×10⁻⁴ mM and the sensor spans its range: **z = 0.028 at rest →
   0.997 during a burst**.

SK2 then does exactly what the review (§V.A.1) describes — it keeps bursts
**short** without abolishing them:

| `gkbar_sk2` | spikes | freq | burst |
|---|---|---|---|
| 0 | 4 | 73 Hz | 41 ms |
| **5×10⁻⁵** (default) | **3** | **63 Hz** | **32 ms** |
| 1×10⁻⁴ | 2 | 40 Hz | 25 ms |
| ≥1×10⁻³ | 1 | — | abolished |

### New: `ihca.mod` — Ca²⁺-dependent I_h

The canonical thalamic spindle terminator (Destexhe et al. 1996), previously
**absent from both cell types**:

```
C  <-> O                 voltage-dependent opening
P0 + 4 Ca <-> P1         Ca2+ binds the regulating factor   (k1ca, k2)
O  + P1   <-> OL         the "locked open" form             (k3p, k4)
I_h = ghbar * (O + 2*OL) * (v - eh)
```

Ca²⁺ from successive rebound bursts locks I_h open → depolarisation → I_T can no
longer de-inactivate → **spindle terminates**; slow unbinding (k2 = 4×10⁻⁴ /ms,
τ ≈ 2.5 s) → **refractory period**. It reads the same private `sk` pool, so the
regulating Ca²⁺ never perturbs the T-current's GHK reversal.

At cell level it works as specified:

| `ghbar` | peak locked-open | V_m drift | bursts fired |
|---|---|---|---|
| 0 | — | +0.3 mV | 2 |
| 1×10⁻⁵ | 0.93 | **+3.4 mV** | 49 |
| **2×10⁻⁵** (default) | 0.90 | **+5.2 mV** | 62 |

I_h acts as both pacemaker (2 → 62 bursts) and depolarising terminator.

### In-network: I_h sets the loop frequency — and reaches the spindle band

`ThalamicNet` now exposes `gsk_re` / `gh_tc` (it previously **hardcoded
`gsk=0.003`**, 60× the calibrated value and inside the burst-abolishing range —
an early "the mechanism over-damps the loop to 1 Hz" reading was that hardcoded
value, not the calibrated defaults).

With the calibrated defaults the network self-oscillates, and **I_h raises the
loop frequency monotonically** — the physiologically expected pacemaker action:

| `gh_tc` | loop frequency | IQR | spectral peak (unconstrained 1–30 Hz) |
|---|---|---|---|
| 0 | 5.1 Hz | 38 ms | 5.0 Hz |
| 5×10⁻⁶ | 6.5 Hz | 5 ms | 6.0 Hz |
| 2×10⁻⁵ | 8.3 Hz | 3 ms | 8.0 Hz |
| 5×10⁻⁵ | 9.3 Hz | 5 ms | 9.0 Hz |
| 2×10⁻⁴ | 7.1 Hz | 193 ms | 3.0 Hz (irregular — broken regime) |
| **4×10⁻⁴** | **13.4 Hz** | **31 ms** | **13.0 Hz** ✅ |

**At `gh_tc = 4×10⁻⁴` the loop runs at 13 Hz — in the 10–15 Hz spindle band for
the first time in this model.** All frequencies above are measured with an
*unconstrained* 1–30 Hz search.

Two honest caveats:

1. **`4×10⁻⁴` is ~20× Destexhe's `ghbar = 2×10⁻⁵`.** Reaching spindle frequency
   by that large a departure from the published conductance needs justification
   before it counts as a result — it may be compensating for something else
   (single compartment, small population, the 1 Hz drive).
2. **The response is non-monotonic** — 2×10⁻⁴ gives a broken, irregular 3 Hz
   regime between the 9 Hz and 13 Hz points. That discontinuity is unexplained
   and should be mapped before trusting the 13 Hz point.

**Event duration is still unsolved**: max 0.18–0.22 s across the whole sweep,
against the review's 0.5–3 s. Frequency and duration are evidently set by
different mechanisms; getting the band right did not lengthen the events.

## MN5 scaling benchmark (job 44671355, 100 ranks) — NEURON scales linearly

Measured on MareNostrum 5, 100 MPI ranks, 1 s of simulation per point:

| cells | wall (s) | ×realtime | s per 1k cells |
|---|---|---|---|
| 100 | 0.3 | 0.28 | 2.78 |
| 500 | 0.4 | 0.45 | 0.90 |
| 2000 | 0.9 | 0.91 | 0.46 |
| **5000** | **2.1** | **2.14** | **0.43** |

**This retires the `N^1.67` figure.** That was measured thread-parallel on one
laptop and does not describe NEURON under MPI. Cost per 1k cells *falls* with N
and flattens (0.46 → 0.43): fixed overhead dominates at small N, and the
marginal cost per cell is near-constant at scale. Between 2000 and 5000 cells
the exponent is **N^0.92** — essentially linear, comparable to NEST's N^1.06.

**Consequence: a 5000-cell HH production run is cheap.** At 2.14× realtime, 200 s
of biological time costs **~7 minutes** of wall clock (plus setup), well inside a
2 h allocation. The earlier advice to hold MN5 time because NEURON-at-scale was
unaffordable was based on the laptop extrapolation and is withdrawn.

⚠️ **One scaling caveat in the model, not the machine.** Convergence follows the
serial model's `k = n/2` rule, so at 5000 cells each neuron integrates ~2500
inputs (~12.5 M synapses) — biologically implausible and O(N²). The `g/k` fan-in
normalisation keeps the *total* conductance constant, so the dynamics should
carry over, but the connectivity rule should become a fixed convergence
(~100–200 sources) before going beyond this size.

The NEST model stays the working reference and the cross-check for the
simulator-agnostic validator; it is no longer the scale-out path, since its
point neurons cannot express the SK2 / Ca-dependent mechanism.

## Cortical column (Aug 2026): L4 / L2/3 / L5 / L6, closing the loop

`cortex_neuron.py` replaces the sinusoidal cortex proxy that
`tc_network_nrn.ThalamicNet._wire_cortex` used with a real conductance-based
column, matching the layer/cell-type structure and connectivity already used
by the NEST reference model (`tc_architecture.py`, `config/network_auditory_hh.yaml`,
Mushtaq et al. 2024 Table 3):

- **`PYCell`** — regular-spiking pyramidal cell, two compartments (soma +
  apical dend). `hh2` fast Na⁺/K⁺ spikes; the dendrite carries `ical` (a
  Reuveni et al. 1993 HVA Ca²⁺ current, opening on every spike upstroke) +
  `cad` (submembrane Ca²⁺ pool) + `sk2` (SK2 Ca²⁺-activated K⁺). **This is the
  same SK2 + Ca²⁺-pool pair already validated in the thalamic RE cell**
  (`mod/sk2.mod`, `mod/cad.mod`), reused here as the cortical
  spike-frequency-adaptation / UP-state-termination current (review sect.
  VI) — the article's SK2 mechanism, generalised from TRN to cortex.
- **`FSCell`** — PV⁺ fast-spiking basket interneuron, `hh2` only (high
  `gkbar`, minimal adaptation): recruits fast feedforward inhibition onto
  pyramidal somata (review sect. VI: TC spindle bursts recruit PV⁺
  interneurons via GluA2-lacking AMPARs).
- **`CorticalColumn`** — L4 / L2/3 / L5 / L6, each an E (`PYCell`) + I
  (`FSCell`) pool (population ratios from `config/network_auditory_hh.yaml`,
  scaled ×0.2 for a fast demo). Wiring: L4→L2/3→L5→L6 feedforward, L6↔L5
  recurrent, E→I feedforward and I→E feedback inhibition within each layer.

**SK2 gives the column a genuine bistable UP/DOWN cycle, not runaway gamma.**
A naively-strong recurrent column (`gkbar_sk2` too low, e.g. 0 – 5e-4) locks
into a permanent self-sustaining oscillation once perturbed — a bang-bang
network with no stable quiescent state over the tested (g_ff, g_rec) range.
At **`gkbar_sk2 = 8e-4`** the column is silent at rest, and a phasic volley
produces a genuine propagating UP-state that self-terminates, instead of
runaway gamma.

## L5 intrinsic bursting: the column's OWN slow-oscillation pacemaker (Aug 2026 v2)

The first cortical-column pass (above) needed an **externally injected**
volley to L4 once per SO cycle — the column had no oscillator of its own,
which is not how the slow oscillation actually originates (it is a cortical
network phenomenon, review sect. VI / Steriade). This is now fixed:

- **`PYCellIB`** (`cortex_neuron.py`) — L5's intrinsically-bursting pyramidal
  cell ("TuftIB", `tc_architecture.py`): `PYCell` + **`mod/inap.mod`**
  (persistent Na⁺, Compte et al. 2003 / Bazhenov-Timofeev cortical SO
  models) + a **slowed** `cad`/`sk2` recovery (`taur` 80 ms → 500 ms). I_NaP
  turns a single spike into a self-sustaining depolarising plateau burst;
  slowing the SK2/Ca²⁺ pool's clearance is what turns that from
  within-burst adaptation into a genuine **~1.4 Hz relaxation oscillator**
  (measured: 12 UP states over 8 s, inter-event interval 694 ms). Half of
  L5's E population is `PYCellIB`, the rest stay regular-spiking `PYCell`
  (`CorticalColumn(ib_frac={"L5": 0.5})`).
- **L5 → L6 is now `g_ff × 6`**, not `g_ff × 1`: L5's IB cells fire a brief,
  sparse burst (a handful of cells, a few ms) once per cycle, not the
  sustained volley the other feedforward links carry, so the ordinary g_ff
  left L6 below threshold most cycles. At ×6 every L5 burst reliably
  reaches L6 — the propagation the corticothalamic "L6 kick" depends on.
- **L5 → L2/3 feedback** (apical-tuft projection) was added so the
  superficial layers can entrain to the L5 SO; L6 ↔ L5 recurrent excitation
  was found to **desynchronise** the IB cells into irregular double-bursting
  and is now off by default (`g_rec = 0`).

```
../.venv-neuron/bin/python cortex_neuron.py
  L5 E: 1.5 Hz/cell mean, 12 UP states, inter-UP-state interval ~694 ms
  L6 E: 2.1 Hz/cell mean, 12 UP states, inter-UP-state interval ~694 ms
  (L4/L2/3 need thalamic/sensory drive to fire -- see ctx_thalamus_network.py)
```

### Bug found and fixed: `mod/ical.mod`'s `vtrap` was inverted

The HVA Ca²⁺ current's alpha-rate helper used `x/(1-exp(-x/y))` instead of
the standard `x/(exp(x/y)-1)` (the form `mod/hh2.mod`'s own `vtrap` uses).
For the Reuveni et al. 1993 alpha formula's operating range (`x/y ≈ 11` at
rest) that is a **~5×10⁴× error**, which pinned the Ca²⁺ channel at ~78%
open at -68 mV instead of near-zero — a lone `PYCell` fired **spontaneously
at ~300 Hz with zero external input**, independent of any bursting
mechanism. Confirmed fixed: a lone `PYCell` (RS, no `gnap`) is now silent at
rest with no input, as it should be.

### The corticothalamic loop: reliable, non-pathological, but not yet a sustained 10–15 Hz train

Two more findings, in order of how the tuning converged:

1. **`gh_tc = 4e-4`** (the thalamus-only "13.4 Hz" result, `ThalamicNet`
   above) makes an **isolated** `TCCell` fire spontaneously at ~103 Hz with
   zero synaptic input — a pathological intrinsic pacemaker, not a rebound
   burst. It only worked in `ThalamicNet` because of one hand-tuned
   synchronous `NetStim` clock driving RE and TC together every cycle; a
   real spiking cortex (propagation delay, jitter) does not reproduce that
   fragile balance, and TC free-ran continuously (~40–100 Hz/cell) once
   coupled to the column. **`ctx_thalamus_network.py` now defaults to
   `gh_tc = 0`** (Destexhe et al.'s Ca²⁺-dependent I_h off), which keeps a
   lone `TCCell` silent with zero input under any condition tested.
2. With `gh_tc = 0`, sweeping the reciprocal RE↔TC gains against the real
   L5-paced L6 kick (not a synchronous clock) found a working point at
   **`g_re_tc = 0.015`, `g_tc_re = 0.011`**: each SO cycle's L6 kick drives a
   clean RE burst (T-current rebound, ~20 spikes over ~60 ms) and 1–5
   TC↔RE volleys at ~140 ms spacing (~7 Hz) before falling silent until the
   next L6 kick 694 ms later — a genuine, non-pathological, SO-locked
   thalamic burst event, and the first multi-cycle (as opposed to
   single-shot) oscillatory train this model has produced.

**What's still short of a real spindle:** the review's definition is a
**sustained 10–15 Hz oscillation for 0.5–3 s** with a waxing/waning envelope.
What's measured here is faster within a burst (individual RE spikes reach
150–500 Hz intra-burst) but the RE↔TC **volley-to-volley** rate is ~7 Hz, not
10–15 Hz, and the number of volleys per event is **inconsistent across
cycles** (4–5 cycles on the first event, tapering to 1–2 on later ones) —
the loop is not yet a clean, reliably-sustained resonator. The population
spectrum (thalamic spike-rate FFT) accordingly peaks at the **SO event rate
(~1.4 Hz)**, not 10–15 Hz; that number describes *how often* an event
occurs, not the oscillation frequency *within* an event.

## Heterogeneity (Aug 2026): a real lever, with a working point and a failure mode

Every TC/RE/cortical cell above is parametrically identical — the only
variability came from which random subset of synaptic partners each cell
drew. That is a plausible reason the RE↔TC loop kept mode-locking to exactly
the SO period (1–2 cycles/event) instead of ringing for several: with only
10 TC + 10 RE cells, a handful of wiring motifs are all there is. Tested
this directly rather than guessing: `ctx_thalamus_network.CorticoThalamicNet`
now takes `het` (uniform ± `het` relative jitter on `e_pas` and `gcabar_it`/
`gcabar_it2` for TC/RE, and on `e_pas`/`gnap` for cortical `PYCellIB`) and
`delay_jitter` (± ms jitter on every synaptic delay), both propagated into
`cortex_neuron.CorticalColumn` too.

**Swept `het` ∈ {0, 0.05, 0.1, 0.2, 0.3} × `delay_jitter` ∈ {0, 1, 2 ms} × 5
seeds, 15 s each (75 runs, ~29 min wall, within the 1 h budget):**

| `het` | events with >3 cycles (of 5 seeds) | mean cycles/event | mean within-event Hz |
|---|---|---|---|
| 0.00 | 5/5 | 2.5 | 7.9 |
| **0.05** | **5/5** | **7.5** | 6.3 |
| 0.10 | 5/5 | 3.0 | 5.9 |
| 0.20 | 2–3/5 | 7.5–8.5 | 5.5–7.6 |
| 0.30 | 0–1/5 | ~1 or degenerate | mostly 0 (collapsed) |

`delay_jitter` made little difference on its own at any `het`. **`het = 0.05`
triples the mean cycle count (2.5 → 7.5) while staying stable across all 5
seeds** — genuine multi-cycle burst trains instead of the near-uniform
1–2-cycle events the homogeneous network produced (`out/` scratch figures:
one seed gave 12 discrete events over 15 s with 4–20 cycles each, vs. the
homogeneous network's near-constant 1–2). **`het ≥ 0.2` starts
destabilising the network**: some seeds still work, but others either
degenerate to a single ~15 s non-SO-locked run (the loop stops resetting
between SO cycles) or nearly full silence — a different, worse failure
mode than the homogeneous network's rigid-but-stable 1–2 cycles. `het =
0.05` is now the default; higher values are available but not recommended
without re-checking each seed individually.

**What this does and does not answer.** Heterogeneity is confirmed as a
real, working lever for event *duration* (cycle count) — not a hypothesis,
a measured 3× effect. It does **not** move the within-event frequency
into the 10–15 Hz band (5.5–7.9 Hz throughout the sweep, no trend toward
15 Hz as `het` rises) — frequency looks set by something else (loop
conduction delay + synaptic time constants, or coarse quantization from
`n_re=10` giving very few effective phase-lag choices) that this cheap,
same-N experiment cannot rule out. That is the argument for actually trying
the bio-plausible cell count next: 5k cells gives a near-continuous
distribution of both intrinsic parameters and RE↔TC delays "for free" from
sampling, at a `het` too fine-grained to reach any other way, which is one
of the two things this local sweep could not test (the other being genuine
population-level progressive recruitment during waxing/waning). See
`ctx_thalamus_mpi.py` / `MN5_NEURON.md` for the scale-out path this
motivates.

## MN5 5k-scale result (Aug 2026): heterogeneity on the L5 pacemaker breaks the loop entirely

The first bio-plausible production run (100 ranks, 5031 cells, `SCALE=1.65`,
`TSTOP=200000`, job 44844450 -- see git history for the full exchange) came
back structurally broken, not just "not yet spindle-band": **RE and
L4/L2-3 never fired even once in 200 s**, L6E fired 26 times total (all in
the first 52 ms -- the init transient, never again), and L5's inhibitory
population (`l5i`) sat at a flat **12 Hz for the entire 200 s** from the
first 5 s bin onward -- not a transient runaway, a stable wrong attractor.

**Root cause, confirmed by a direct A/B test** (`ctx_thalamus_mpi.py`,
scale=0.3, 159 L5E cells, 8 s): with `het` jittering `PYCellIB`'s burst-timing
parameters (`e_pas`, `gnap`), L6E and RE both stay at **0 active cells**;
turning that jitter off (`het=0.0` everywhere) restores **126/126 L6E and
16/16 RE active**. The Aug 2026 heterogeneity study above only tested this at
n≈10 thalamic cells and a ~5-cell L5 IB subpopulation -- too small to expose
the failure mode. At bio-plausible scale (400+ independent IB pacemakers),
jittering each cell's burst phase/threshold is exactly a weakly-coupled-
oscillator desynchronisation problem: without shared timing, the population
drifts out of phase and never produces a coincident burst again after the
initial (artificial, initialisation-driven) transient. No coincident L5
burst means L5→L6 never crosses threshold, which cascades: no L6 output →
no corticothalamic drive → RE never triggered → TC never rebounds → L4/L2-3
(which have no intrinsic driver of their own) never wake up. Meanwhile L5's
interneurons, receiving a smeared-out, now-continuous-rather-than-phasic
barrage from 100 (`conv`) independently-timed sources, settle into tonic
firing instead of a burst terminator.

**Fix, verified by the same A/B test with the fix applied:** keep
heterogeneity everywhere it was shown to help (TC, RE, cortical RS/FS cells)
but make `PYCellIB`'s `e_pas`/`gnap` homogeneous again -- the SO-generating
pacemaker population needs to stay in phase; nothing downstream of it needs
it to be individually jittered. With the fix, the same 159-L5E-cell network
at `het=0.05` reaches **100% active L6E (126/126), RE (16/16), and TC
(63/63)** -- matching or beating the `het=0` control. Both
`cortex_neuron.CorticalColumn` and `ctx_thalamus_mpi.ParallelCorticoThalamicNet`
are fixed; the fix has NOT yet been re-run at the full 5031-cell / 200 s MN5
production scale (that is the next thing to do once this is on `main`).

**Lesson for future scale-out work in this model:** "does bio-plausible N
fix an open dynamical question" is not a safe default hope -- a mechanism
that depends on population-level synchrony (here, the SO pacemaker) can
break, not improve, when naively scaled with per-cell heterogeneity, even
though heterogeneity was independently validated as beneficial for a
*different* part of the same network (the RE↔TC loop). Test each population
separately before assuming a global heterogeneity/scale knob helps
uniformly.

### Known limitations of this cortical column

- Population sizes are scaled ×0.2 from `config/network_auditory_hh.yaml`
  for demo speed; not yet validated at full scale.
- No matrix (diffuse, TC→L1) pathway — only the core (TC→L4, focal) pathway
  is wired, appropriate for a first-order auditory column (review sect. V.B
  / `docs/spindle_review_mapping.md` §3.5) but meaning only fast-spindle-type
  dynamics are representable here.
- No dendritic PV→SST disinhibition (review sect. VI) — `FSCell` inhibits
  `PYCell` somata directly; SST⁺ interneurons and dendritic Ca²⁺
  disinhibition are out of scope for this pass (`docs/spindle_review_mapping.md`
  §3.6 notes the same limitation for the NEST model).
- The corticothalamic/thalamocortical conductances (`g_tc_l4`, `g_l6_tc`,
  `g_l6_re`) and the column's own recurrent gains (`g_ff`, `g_rec`, `gkbar_sk2`)
  were hand-tuned to a stable, non-pathological operating point, not swept
  systematically — treat them as a starting point for further calibration,
  not a validated result.
- The RE↔TC loop produces a variable number of oscillatory cycles per SO
  event (4–5 on the first event, 1–2 on later ones) at ~7 Hz, not a
  consistent 10–15 Hz train lasting 0.5–3 s — see "not yet a sustained
  10–15 Hz train" above. This is the main open item before calling the
  model's output a spindle in the review's strict sense.
