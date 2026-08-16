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
| `tc_neuron.py` | `TCCell` (relay), `RECell` (reticular), `gap_junction()` + demos |
| `tc_network_nrn.py` | `ThalamicNet` — the TC↔RE loop network; emits the `tc_validate` contract |
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
