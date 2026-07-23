# Mapping Fernandez & Lüthi (2020) onto the auditory thalamo-cortical model

Reference: **Fernandez LMJ, Lüthi A. "Sleep Spindles: Mechanisms and Functions."**
*Physiol Rev* 100: 805–868, 2020.
[doi:10.1152/physrev.00042.2018](https://doi.org/10.1152/physrev.00042.2018)

This note summarises the review and records, mechanism by mechanism, what our
model already captures, what is missing, and the concrete consequences for
future work. It also proposes designs for two specific extensions:
**SK2 channels in TRN** (§4) and an **external modulatory spindle trigger** (§5).

---

## 1. Summary of the review

The review argues for an **updated view**: sleep spindles are not merely a
sinusoidal EEG marker but a **predominantly local thalamocortical event** that

1. times NREMS **continuity vs. fragility**,
2. reshapes the cortical **excitation–inhibition balance**,
3. enables **active memory consolidation** (reactivation/routing between
   hippocampus and cortex),
4. serves as a **biomarker** for cognitive impairment and schizophrenia.

### Definitions and measurement (§III)

- Spindles: trains of **10–15 Hz** cycles (AASM: 11–16 Hz, most 12–14 Hz),
  duration ≥0.5 s, waxing/waning envelope; hallmark of NREMS stage N2.
- Quantified either as **sigma power** (7–16 Hz, typically 10–15 Hz) or as
  **discrete events** (fixed thresholding, or Morlet wavelet on a 9–16 Hz band).
- Key metrics: **density** (number per unit time) and **intra-spindle
  frequency** (cycle number / duration).
- **Fast** spindles (~13–15 Hz, posterior/sensory) vs **slow** (~10–12 Hz,
  frontal/anterior); the border is conventionally ~13 Hz.

### Circuit mechanism (§V)

The **thalamic reticular nucleus (TRN, our nRT/RE)** is the spindle pacemaker,
operating in a reciprocal loop with **thalamocortical relay (TC) cells**.

| phase | mechanism |
|-------|-----------|
| **Initiation** | TRN fires **low-threshold Ca²⁺ bursts** (Ca_v3.3 / T-type), possible **only when hyperpolarised below ~−55 mV**. **SK2** (Ca²⁺-activated K⁺) channels produce the burst after-hyperpolarisation, keeping bursts short and repeatable. Triggered by **layer-6 corticothalamic** afferents (fast GluA4 AMPARs) and by TC→TRN excitation. |
| **Inhibition** | TRN→TC via **GABA_A (synaptic α1γ2 + extrasynaptic α4/δ) *and* GABA_B**. Burst-driven IPSPs are disproportionately large and slow (transmitter spillover recruits extrasynaptic receptors). A **tonic extrasynaptic δ-GABA_A** current holds TC hyperpolarised and *facilitates* rebound. |
| **Rebound** | TC rebound bursts via **Ca_v3.1** (primary dendrites), only below ~**−65 mV**, and typically only after **several IPSPs** — the origin of the spindle's *waxing* build-up. |
| **Synchronisation** | *Promoted* by topographic L6 corticothalamic **"slabs"**, **gap junctions** (connexin-36) between TRN cells, and open-loop TC–TRN lateral excitation. *Constrained* by **TRN–TRN lateral inhibition** (GABA_A/GABA_B), which **antisynchronises**. Net: synchrony is **modest**, spindles are **mostly local**. |
| **Termination / refractory (5–10 s)** | TC **I_h (HCN)** after-depolarisation suppresses further rebound; progressive TRN hyperpolarisation via Ca²⁺- and Na⁺-dependent K⁺ currents; **noradrenaline** from locus coeruleus rises toward spindle end. |

**Heterogeneity.** First-order **sensory TRN sectors (explicitly including
auditory)** burst vigorously → strong, focal spindles; higher-order sectors are
weaker. **PV⁺-TRN** cells are focal/first-order and burst strongly;
**SST⁺-TRN** cells project diffusely (candidate "global" spindles).

### Cortical state (§VI)

TC spindle bursts recruit **fast feedforward PV⁺ interneurons** (GluA2-lacking
AMPARs, high release probability), which **inhibit pyramidal somata** while
**disinhibiting dendrites** (PV→SST), raising **dendritic Ca²⁺** — a plasticity
gate. TC **core** projections target L4 focally (first-order); **matrix**
projections target L1 diffusely. Spindles nest on the **slow-oscillation (SO)
up-state** (50–70% phase-locked) and couple to **hippocampal ripples**.

### Temporal organisation (§IV)

Spindles cluster on a **~0.02 Hz (~50 s) infraslow rhythm**, alternating
**continuity** (spindle-rich, arousal-protected) and **fragility** (spindle-poor,
arousable) periods, co-modulated with heart rate, pupil diameter and brain
temperature. Homeostatically, **sigma power is *down*-regulated** by sleep
deprivation (opposite to slow-wave activity) and increased after napping.

---

## 2. What our model already captures

| review mechanism | our implementation | status |
|---|---|---|
| TRN burst pacemaker, hyperpolarisation-gated I_T | `ht_neuron` RE with `g_peak_T`, hyperpolarised via `g_KL_thalamus` | ✅ aligned |
| TRN→TC **GABA_A + GABA_B** | `connect_inh(..., gaba_b=...)` on nRT→TCR | ✅ aligned |
| TC rebound needs hyperpolarisation + several IPSPs | matches our calibration (TC needs hyperpolarised baseline; rebound weak) | ✅ aligned |
| TRN–TRN lateral inhibition (antisynchronising) | nRT→nRT self-inhibition | ✅ aligned |
| L6 corticothalamic feedback initiates spindles | L6→TCR / L6→nRT (2:1, Mushtaq Table 3) | ✅ aligned |
| SO up-state gating of spindles | light 1 Hz scaffold + CT feedback | ✅ aligned |
| First-order **auditory** sensory loop, core (→L4) projection | MGB→L4 thalamocortical input | ✅ aligned |
| Feedforward interneuron recruitment | TC→L4_I projection | ⚠️ partial (see §3.6) |

---

## 3. Key consequences (prioritised backlog)

### 3.1 Enforce a 5–10 s spindle refractory period — *highest realism gap*
The review sets the refractory at **5–10 s**; we currently emit roughly **one
spindle per ~1 s SO cycle**, which is far too frequent. **Consequence:**
strengthen the TC **I_h** after-depolarisation (`g_peak_h`) and/or add a slow
activity-dependent K⁺ current so that *not every* UP state yields a spindle.
Validation metric: inter-spindle interval distribution should peak at 5–10 s.

### 3.2 Add the ~0.02 Hz infraslow clustering (continuity/fragility)
Absent from our model. **Consequence:** modulate thalamic excitability on a
**~50 s** timescale (slow drift of `g_KL` or a tonic bias) so spindles group
into clusters separated by spindle-free periods. This directly implements the
review's function #1 and yields a new validation metric: **spindle density
itself oscillating at ~0.02 Hz**. See §5 — the modulatory-signal design provides
this naturally.

### 3.3 Make thalamic connectivity topographic, not all-to-all
Our intra-thalamic loop uses `prob=1.0` (every RE inhibits every TC), which
**over-synchronises** relative to the review's *modest, local* synchrony.
**Consequence:** replace all-to-all with **topographic/local** connectivity
("slabs"), keeping RE→RE lateral inhibition as the antisynchroniser. Optionally
add **gap junctions between RE cells** (`ht_neuron` supports NEST gap junctions)
for within-sector synchrony. Expected outcome: **local** spindles, more
biological amplitude variability.

### 3.4 SK2 (Ca²⁺-activated K⁺) in TRN
Currently thalamic `g_peak_KNa = 0`, so TRN bursts have no proper
after-hyperpolarisation. **Consequence:** add an SK2-like current to shape burst
duration → affects spindle amplitude and regularity. Full design in **§4**.

### 3.5 Fast vs slow spindles = core (L4) vs matrix (L1)
We only model the core pathway. **Consequence:** to represent both spindle
types, add a diffuse matrix-like TC→L1/upper-layer projection. For the
**auditory** (first-order) column the core pathway should dominate, so our
target band is the **fast** spindle range (~13–15 Hz).

### 3.6 Feedforward PV inhibition and the cortical spindle state
**Consequence:** wire **TC→PV (Basket)** explicitly as the fast feedforward path
and validate that during spindles **pyramidal soma firing is suppressed /
phase-locked while interneuron firing increases**.
**Limitation:** point neurons cannot represent **dendritic disinhibition,
dendritic Ca²⁺, or the CSD depth profile** — explicitly out of scope.

### 3.7 Align measurement with the review
**Consequence:** narrow spindle detection to **10–15 Hz** (fast 13–15 / slow
10–12), and report **spindle density (#/min)** and **intra-spindle frequency**
alongside the existing waxing/waning envelope in the decomposition figure.

### 3.8 Honest scope limits
Point neurons (no dendrites, no CSD, no dendritic Ca²⁺ plasticity); generic
`ht_neuron` I_T/I_h rather than specific Ca_v3.3 / Ca_v3.1 / SK2 / HCN subunits;
no noradrenergic termination; no gap junctions; no infraslow clustering (yet).

---

## 4. Proposed design — SK2 channels in TRN under HH formalism

**Biological role (review §V.A.1).** TRN dendrites express Ca_v3.3 T-type
channels; Ca²⁺ entering through them activates **SK2** channels, producing a
**burst after-hyperpolarisation (AHP)**. This keeps each burst **short** and
allows bursts to **repeat reliably** — directly shaping sleep-spindle amplitude
and regularity. Reduced SK2 function is implicated in attentional disorders and
Dravet syndrome.

### 4.1 Equations (Hodgkin–Huxley form)

SK2 is **voltage-independent** and gated purely by intracellular Ca²⁺, so the
activation variable depends on $[\mathrm{Ca}^{2+}]_i$ rather than $V$:

$$
I_{\mathrm{SK}} = g_{\mathrm{SK}}\; z \;\bigl(V - E_K\bigr)
$$

with a **Hill-type** Ca²⁺ activation gate

$$
z_\infty\bigl([\mathrm{Ca}^{2+}]_i\bigr) =
\frac{[\mathrm{Ca}^{2+}]_i^{\,q}}{[\mathrm{Ca}^{2+}]_i^{\,q} + K_d^{\,q}},
\qquad
\frac{dz}{dt} = \frac{z_\infty - z}{\tau_z}
$$

Typical values: Hill coefficient $q \approx 4$, $K_d \approx 0.3$–$0.7\ \mu M$,
$\tau_z \approx 5$–$15$ ms (fast activation, slower deactivation), $E_K = -90$
to $-95$ mV.

The gate needs a **fast Ca²⁺ pool** driven by the T-type current:

$$
\frac{d[\mathrm{Ca}^{2+}]_i}{dt}
  = -\beta_{\mathrm{Ca}}\, I_T \;-\; \frac{[\mathrm{Ca}^{2+}]_i - [\mathrm{Ca}^{2+}]_{\mathrm{rest}}}{\tau_{\mathrm{Ca}}}
$$

with $\tau_{\mathrm{Ca}} \approx 10$–$50$ ms for the SK-relevant submembrane pool.

### 4.2 Implementation options in NEST

| option | effort | fidelity | notes |
|---|---|---|---|
| **(a) I_KNa as an SK2 surrogate** | trivial (set `g_peak_KNa > 0` on RE) | approximate | `ht_neuron` already has a **Na⁺**-dependent K⁺ current. It is *gated by Na⁺, not Ca²⁺*, but is likewise **spike/burst-activated** and produces an AHP, so it reproduces the *functional* role (burst shortening, repeatable bursting). **Recommended first step** — testable today, no new model. |
| **(b) NESTML custom model** | moderate | high | Generate `ht_neuron_sk` with NESTML, adding $I_{\mathrm{SK}}$ and a fast Ca²⁺ pool. Correct Ca²⁺ gating and a genuine SK2 parameterisation. **Recommended proper solution.** |
| **(c) C++ NEST module** | heavy | high | Only if NESTML proves limiting. |

**Important caveat for (b):** `ht_neuron` *already* maintains a Ca²⁺ variable
(`beta_Ca`, `tau_Ca`) fed by I_T — but its default `tau_Ca = 10000 ms` is tuned
for the slow **Ca²⁺–cAMP–HCN** mechanism behind spindle termination. SK2 needs a
**much faster** transient (tens of ms), so a custom model must maintain a
**second, fast Ca²⁺ pool** for SK rather than reusing the existing slow one.

### 4.3 NESTML sketch

```
neuron ht_neuron_sk:
  state:
    z real = 0.0                  # SK2 activation gate
    Ca_fast uM = 50.0 nM          # fast submembrane Ca pool (SK-relevant)
  parameters:
    g_SK   nS  = 2.0
    E_K    mV  = -90.0 mV
    K_d    uM  = 0.5 uM
    q      real = 4.0
    tau_z  ms  = 10.0 ms
    tau_Ca_fast ms = 20.0 ms
    beta_Ca_fast real = 5e-3
  equations:
    # fast Ca pool driven by the T-type current
    Ca_fast' = -beta_Ca_fast * I_T - (Ca_fast - Ca_rest) / tau_Ca_fast
    # Hill-type Ca gating (voltage independent)
    inline z_inf real = Ca_fast**q / (Ca_fast**q + K_d**q)
    z' = (z_inf - z) / tau_z
    inline I_SK pA = g_SK * z * (V_m - E_K)
    # ... add I_SK to the membrane current sum
```

### 4.4 Implemented (stage 1) — measured result

Stage (a) is **implemented**: `HHParams.g_peak_KNa_RE` (default `0.5`) switches
`ht_neuron`'s Na⁺-dependent K⁺ current on **for the RE/nRT population only**.
Sweeping it reproduces the review's prediction that the burst AHP keeps TRN
bursts **short and repeatable** (20 s runs, local config):

| `g_peak_KNa_RE` | RE spikes/burst | number of bursts | spindle-envelope CV |
|---|---|---|---|
| 0.0 | 5.97 | 144 | 0.876 |
| 0.5 | 5.08 | 151 | 0.859 |
| 1.0 | 4.93 | 151 | 0.878 |

Bursts shorten by ~15% and become slightly **more frequent** — i.e. shorter and
more repeatable, as the review describes. **Honest caveat:** the predicted
improvement in *spindle regularity* is **not** demonstrated (envelope CV is
essentially unchanged); that may require the true Ca²⁺-gated SK2 of option (b),
or regularity may be dominated by other factors in our circuit.

### 4.5 Where it plugs into our model

Apply **only to the reticular (RE/nRT) population** — per the review, SK2
shapes *TRN* burst discharge (`_ht_neuron_params(role="RE")` in
[`tc_network.py`](../tc_network.py)). Expected effects, and how to validate:

- **Shorter, more stereotyped RE bursts** → measure spikes-per-burst.
- **More regular spindle envelope** → lower variance of the Hilbert envelope.
- **Modulated spindle amplitude/duration** → $g_{\mathrm{SK}}$ becomes a
  principled knob for spindle regularity (and a disease model: reduced SK2 →
  degraded spindles, cf. review §IX/X).

---

## 5. Proposed design — external modulatory signal to trigger spindles

**The distinction that matters.** Triggering a spindle externally is
biologically correct — the review is explicit that spindles are *initiated* by
**layer-6 corticothalamic** input and *gated* by **neuromodulation**
(noradrenaline terminates them; melatonin promotes TRN bursting; the ~0.02 Hz
infraslow rhythm groups them). What must **not** be imposed externally is the
**intra-spindle frequency** — that has to stay emergent from the RE↔TC loop, or
we are back to the old imposed-oscillator model. So:

> **Externally trigger *when* a spindle happens; let the loop decide *at what
> frequency* and *for how long*.**

Our current HH model already follows the good half of this rule (light 1 Hz
scaffold gates *when*; the loop sets the frequency). The proposal below makes
the trigger explicit and adds the missing slow gating.

### 5.1 Two-component modulatory signal

**(A) Phasic trigger — the corticothalamic "kick" (initiation).**
A brief depolarising pulse (or a small spike volley) delivered to **RE/nRT** at
the SO up-state onset, standing in for the L6 corticothalamic burst.
Implementation: `dc_generator`/`step_current_generator` or a `spike_generator`
→ AMPA onto RE, phase-locked to the SO. Duration ~20–50 ms; the spindle that
follows is generated by the loop, not by the pulse.

**(B) Tonic/slow modulator — sleep depth and clustering (permission).**
A slow signal that sets **whether a triggered spindle succeeds**, mapping onto
Hill & Tononi's explicit neuromodulation knob **`g_KL`** (K-leak) — larger
`g_KL` = deeper hyperpolarisation = more burst-permissive. Give it a **~0.02 Hz
(~50 s)** component and the model reproduces **continuity vs fragility**
periods (§3.2) for free. A noradrenaline-like *rise* can additionally be used to
**terminate** spindles and impose the **5–10 s refractory** (§3.1).

### 5.2 Why this is the right architecture

- It implements **three backlog items at once**: initiation (§5.1A), infraslow
  clustering (§3.2), and refractoriness/termination (§3.1).
- It keeps the scientifically important part — **intra-spindle frequency,
  waxing/waning, and duration** — **emergent**.
- It is directly interpretable: component (A) is corticothalamic drive,
  component (B) is neuromodulatory tone. Both are named mechanisms in the review.

### 5.3 Suggested parameterisation

| signal | target | form | rationale |
|---|---|---|---|
| phasic trigger | RE/nRT | pulse at SO up-state onset, 20–50 ms | L6 CT burst (review §V.B.2) |
| slow permission | TC + RE | `g_KL` (or tonic bias) modulated at **0.02 Hz** | continuity/fragility (§IV.E) |
| termination | TC + RE | transient depolarising/NA-like rise after spindle onset | LC-noradrenaline termination (§V.D.3) |

### 5.4 Implemented — and what the calibration taught us

Both components are **implemented** for the HH model:

- **(A) phasic trigger** — `SleepParams.spindle_trigger` (CLI `--spindle-trigger`)
  attaches a `step_current_generator` delivering a `trigger_amp` (30 pA),
  `trigger_dur` (30 ms) kick to **nRT** once per SO cycle.
- **(B) slow permission** — `HHParams.infraslow_freq` / `infraslow_amp`, a
  0.02 Hz `ac_generator` on the thalamus.

Two findings from calibrating this, both worth recording:

1. **The infraslow swing must be large (~40 pA), and its effect is
   non-monotonic.** Hyperpolarising the thalamus lowers TC firing but
   *simultaneously* makes RE **more** burst-prone (I_T de-inactivates), so the
   two effects partly cancel. Measured sigma-envelope modulation ratio at the
   infraslow frequency: **3.1** (no modulation) → **2.6** at 22 pA (no effect)
   → **8.8** at 40 pA. Hence the 40 pA default.
2. **Discrete-event counts cannot validate clustering.** Our `detect_sp_sw`
   uses a *percentile* threshold, so it self-normalises and reports a similar
   event count however strongly spindles are modulated. The review's own
   formulation — "**the ~0.02 Hz oscillation of sigma power**" — is the correct
   measure, so the runner now reports a **`sigma infraslow mod`** metric: the
   spectral peak of the spindle-band envelope inside the infraslow band, with a
   power ratio. This is the metric to trust for §3.2.

**Validation run** (local config, `ht_neuron`, `--spindle-trigger`, 200 s ≈ 4
infraslow cycles):

| metric | result | target |
|---|---|---|
| slow-wave peak | **1.00 Hz** | ~1 Hz ✅ |
| spindle peak | **11.0 Hz** | 10–15 Hz ✅ (emergent, not imposed) |
| **sigma infraslow modulation** | **peak 0.020 Hz, ratio 14.6** | ~0.02 Hz clustering ✅ |
| inter-spindle interval | 0.98 s | 5–10 s ❌ (see below) |

The spindle frequency stays in the physiological band **while the trigger runs
at 1 Hz** — confirming the intra-spindle rhythm is set by the RE↔TC loop, not by
the modulatory signal, which was the whole design constraint.

### 5.4.1 The 5–10 s refractory: two mechanisms tried, both insufficient

**Status: not achieved.** Two principled attempts, both measured with an
*absolute* spindle metric (contiguous runs of TC firing — see the note below on
why the percentile detector cannot be used here):

**Attempt 1 — deepen the TC I_h after-depolarisation** (the review's own
termination mechanism: Ca²⁺ via Ca_v3.1 → cAMP → HCN). `HHParams.beta_Ca_TC`
and `g_peak_h_thalamus` are now exposed as levers. Sweeping them *does* lengthen
the inter-spindle interval, but at an unacceptable cost:

| `g_peak_h` | `beta_Ca_TC` | episodes/40 s | median ISI | duration | spindle freq |
|---|---|---|---|---|---|
| 2 | 0.001 | 9 | 2.98 s | — | 11–15 Hz ✅ |
| 3 | 0.003 | 7 | **4.99 s** | 5.2 s ❌ | <9 Hz ❌ |
| 6 | 0.005 | 3 | 17.4 s ❌ | 5.2 s ❌ | <9 Hz ❌ |
| 12 | 0.05 | 2 | — | tonic ❌ | destroyed |

Pushing I_h far enough to space spindles out **drags the intra-spindle frequency
below 9 Hz** (out of the physiological band) and **stretches episodes past 3 s**
(the review says 0.5–3 s); pushed further, TC leaves burst mode entirely and
fires tonically, abolishing spindles. There is a genuine trade-off here in
`ht_neuron`: the same I_h that terminates a spindle also depolarises the relay.

**Attempt 2 — space the external trigger** (`trigger_refractory_ms`), motivated
by the fact that human spindle density is only ~2–8/min, so the corticothalamic
trigger cannot be firing on every SO cycle. **Measured: no effect** — with the
trigger at 7 s the model still produced **40 episodes in 60 s (ISI 1.03 s)**.
The reason is instructive: the **1 Hz SO drive plus L6 corticothalamic feedback
re-activate the thalamus on every UP state regardless of the trigger**, so the
trigger is additive, not the sole initiator.

**What this implies for the next attempt.** A refractory requires making the
thalamus *non-permissive between spindles*, not just triggering it less often.
Candidates: (i) a much deeper/slower permission signal so most UP states find
the thalamus out of burst range; (ii) reducing the tonic SO drive to the
thalamus so *only* triggered cycles spindle; (iii) the true Ca²⁺-gated SK2 of
§4 option (b), whose AHP acts on RE without depolarising TC the way I_h does.

**Measurement caveat (important).** The percentile-based `detect_sp_sw` cannot
be used to test a refractory, for the same reason it could not test clustering:
it self-normalises. It reported ISI ≈ 0.98 s under settings where the absolute
TC-episode measure reported ≈ 3 s. Use the absolute measure
(`tc_episodes()` in [`tc_spindle_figures.py`](../tc_spindle_figures.py)).

**Previously reported as still open:** the **5–10 s refractory (§3.1)** — spindles
still recur about once per ~1 s SO cycle (median inter-spindle interval
≈ 0.95 s). Deepening the TC I_h after-depolarisation (`g_peak_h`, and the slow
Ca²⁺–cAMP pathway that `ht_neuron` already models with `tau_Ca = 10000 ms`)
is the natural next lever, possibly combined with triggering on only a subset of
SO cycles.

### 5.5 Validation metrics

1. **Intra-spindle frequency** stays in 10–15 Hz and is *independent* of the
   trigger rate (proves the frequency is emergent, not imposed).
2. **Inter-spindle interval** shows a 5–10 s refractory.
3. **Spindle density** oscillates at ~0.02 Hz (clustering).
4. **SO-phase locking**: 50–70% of spindles on the SO up-state.

---

## 5.4.2 Gated mode: third refractory attempt, also unsuccessful

Rationale: attempts 1–2 failed because the thalamus is re-activated on **every**
SO cycle, so spindle density is pinned near the SO rate. "Gated mode"
(`SleepParams.trigger_refractory_ms > 0`) therefore holds the relay
**sub-threshold between spindles** (`HHParams.tc_slow_offset_gated`), makes the
infraslow signal **purely suppressive**, and lets the trigger open a permissive
**window** (`trigger_window_ms`, `trigger_tc_amp`) that sets spindle duration.

**Result: does not work, and is off by default.** A first pass
(`offset 18 + SO amp 12 = 30 pA`) produced *no gating at all* — 30 pA is exactly
the level at which TC fires under reticular inhibition, so every SO peak still
re-fired the relay. Closing the gate properly (`offset 6`, stronger window)
instead **suppressed the thalamus**: MGB 4.0 → 1.7 Hz, nRT 10 → 2.7 Hz, and the
loop fell to ~9 Hz, failing the runner's own spindle check — while spindle
density did **not** fall. The parameters are retained as an explicit opt-in
(`trigger_refractory_ms = 0` by default, restoring the validated behaviour).

## 5.5 Verdict: do the simulated events qualify as spindles?

**No — not by the review's own definition.** Measured with the review's
detection method (fixed threshold on the 10–15 Hz envelope, events ≥ 0.5 s;
`detect_spindles()` in [`tc_validate.py`](../tc_validate.py)):

| criterion | paper | measured | |
|---|---|---|---|
| intra-spindle frequency | 10–15 Hz | 15.1 Hz | ✗ fast |
| **duration** | **0.5–3 s** | **max 0.45 s, median ~0.2 s** | **✗** |
| density (events ≥0.5 s) | 2–8 /min | **0.0 /min** | ✗ |
| inter-spindle interval | 5–10 s | ~1 s | ✗ |
| SO coupling | 50–70% | 35–86% (unstable) | ✗ |
| infraslow clustering | ~0.02 Hz | 0.020 Hz, ratio 15 | ✓ |
| **RE spikes per burst** | **2 to >10 (>100 Hz)** | **1.09–1.14** | **✗** (see §5.5.1) |
| **TC spikes per burst** | **2–6 (rebound burst)** | **1.21–1.31** | **✗** (see §5.5.1) |
| RE V_m < −55 mV | required | 72–94% of time | ✓ |
| TC V_m < −65 mV | required | 47–73% of time | ✓ |

> The two "spikes per burst" rows previously read as passes. They were measured
> with a 20 ms intra-burst window, which lumps tonic spikes into pseudo-bursts;
> at the review's own >100 Hz criterion they fail. See §5.5.1 — this is the
> project's key negative finding.

**The finding is robust, not a detection artefact:** across thresholds
(1.0/1.5/2.0 SD) and in both gated and ungated modes, the **longest**
sigma-band event is **0.45 s** and **none reaches 0.5 s**.

**Interpretation.** The model reproduces the **cellular and circuit mechanism**
faithfully — correct membrane-potential operating ranges, realistic TRN burst
structure, a working RE↔TC loop, and correct infraslow organisation. What it
does **not** yet reproduce is the **event**: a spindle is a *train* of 7–45
cycles (0.5–3 s at 10–15 Hz), whereas our loop rings for only ~2–3 cycles
(~0.2 s) before dying. It is best described as **sigma-band ringing, not sleep
spindles**.

### 5.5.1 Root cause identified: the cells fire tonically, not in bursts

Follow-up measurement (see §5.6) **corrects one of the passes above** and
identifies the real blocker.

At the review's own burst criterion — a burst is 2 to >10 spikes at up to
several hundred Hz, i.e. **intra-burst ISI < 10 ms** — the model produces:

| | measured | review |
|---|---|---|
| TC spikes per burst | **1.21** | 2–6 (rebound burst) |
| RE spikes per burst | **1.09** | 2 to >10 |

That is **single spikes: the thalamus is in TONIC firing mode, not BURST
mode.** The earlier "RE spikes per burst 3.6–5.2 → PASS" used a lenient 20 ms
window, which lumps unrelated tonic spikes into pseudo-bursts; at the correct
>100 Hz criterion it fails. The validator now uses the 10 ms criterion.

TC sits at mean −62.7 mV (p95 −54.1), i.e. hovering just under threshold and
emitting isolated spikes, rather than resting hyperpolarised and firing
**post-inhibitory rebound bursts** — the actual spindle mechanism.

This traces back to the very first `ht_neuron` calibration in this project:
hyperpolarising a relay cell and releasing it produced a rebound that peaked
near **−56 mV, below the −51 mV threshold** — a sub-threshold rebound, and
raising `g_peak_T` (1 → 3 → 6) barely changed it. The workaround at the time was
to depolarise TC until it fired, which yields *tonic firing gated by rhythmic
inhibition* — superficially spindle-like in the LFP, but not the burst mechanism.

**Consequence: the blocker is upstream of connectivity.** No wiring change can
sustain a spindle if the relay never emits a rebound burst — which is why the
topographic/open-loop experiment below had no effect. Fixing this requires a
neuron model whose T-current genuinely produces rebound bursts (a NESTML model
with an explicit Ca_v3.1 low-threshold current, cf. §4 option b), not further
tuning of `ht_neuron`.

### 5.5.2 Topographic wiring and open-loop recruitment: implemented, no effect

Implemented (`HHParams.topographic`, opt-in): ring-topographic RE→TC ("slabs"),
local RE→RE lateral inhibition, and the review's **open-loop** TC→RE projection
(`tc_re_offset`), where a relay cell excites reticular cells *adjacent* to those
that inhibited it — the mechanism the review credits with laterally propagating
and sustaining spindle rhythms.

**Measured: no effect on spindle duration** (40 s runs, 40 TC / 40 RE):

| wiring | events | median dur | max dur | events ≥0.5 s |
|---|---|---|---|---|
| all-to-all (baseline) | 22 | 0.18 s | 0.24 s | 0 |
| topographic, offset 0 | 18 | 0.18 s | 0.24 s | 0 |
| topographic, offset 2 | 17 | 0.21 s | 0.27 s | 0 |
| topographic, offset 4 | 18 | 0.17 s | 0.24 s | 0 |

Consistent with §5.5.1: recruitment cannot sustain an oscillation whose
constituent cells never burst. Kept as an opt-in (`topographic = False`) for
re-testing once a burst-capable relay model exists.

**Gap junctions — not implementable.** NEST's `gap_junction` synapse requires
waveform relaxation; `ht_neuron` reports `node_uses_wfr = False`, and the only
gap-junction-capable model available here (`hh_psc_alpha_gap`) lacks I_T/I_h.
Verified by direct test (`IllegalConnection`).

**Superseded root-cause hypothesis.** The loop is under-damped for too short a
time —
it does not reverberate long enough. Per the review's synchronisation section
(V.C), sustained spindles depend on mechanisms we have not implemented:
**gap junctions between TRN cells**, **topographic corticothalamic "slabs"**,
and **open-loop TC→TRN lateral excitation** that recruits progressively more
TRN cells over the course of a spindle. Our all-to-all thalamic wiring (§3.3)
has no such recruitment dynamics, so there is nothing to sustain the
oscillation. Implementing §3.3 is therefore the highest-value next step — it is
likely the prerequisite for spindle *duration*, and hence for density and the
refractory period too.

## 5.6 Resolution: AdEx thalamus produces real spindles (6/8 criteria)

The blocker in §5.5.1 was that `ht_neuron`'s I_T never fires a rebound burst.
The fix is a relay model that does. **NESTML is not installed here**, but a
cheaper option was already in the model list: **`aeif_cond_exp` (AdEx)**.

**Mechanism.** AdEx reproduces the post-inhibitory rebound through its
subthreshold adaptation. While the cell is hyperpolarised, $a\,(V-E_L)$ drives
$w$ negative; on release $-w$ is **depolarising**, and if it carries $V$ to
$V_{th}$ the exponential term regenerates and fires a high-frequency burst,
which $b$ then terminates.

**Single-cell calibration.** A scan over $V_{th}$, $g_L$, $\Delta_T$, $a$, $b$,
$\tau_w$, $V_{reset}$ found 48 rebound-bursting parameter sets. The chosen one
(`AdExParams`, thalamic block) gives **6 spikes at ~298 Hz** after release from
hyperpolarisation, silent at rest and silent afterwards — matching the review's
2–6 spike rebound burst. (Note the first attempt failed: with the default
$V_{th}=-50$ mV the rebound peaked at −57 mV, the *same* 7 mV shortfall as
`ht_neuron`; lowering $V_{th}$ and $g_L$ is what closes it.)

**Network result** (40 TC / 40 RE, 60 s), versus the same test on `ht_neuron`:

| criterion | paper | ht_neuron | **AdEx** |
|---|---|---|---|
| TC spikes per burst | 2–6 (>100 Hz) | 1.21 ✗ | **4.69 ✓** |
| RE spikes per burst | 2 to >10 | 1.09 ✗ | **3.92 ✓** |
| duration | 0.5–3 s | 0.2 s (max 0.45) ✗ | **0.51 s ✓** |
| intra-spindle frequency | 10–15 Hz | 15.1 ✗ | **13.7 ✓** |
| RE V_m < −55 mV | required | ✓ | ✓ |
| TC V_m < −65 mV | required | ✓ | ✓ |
| density | 2–8 /min | 0.0 ✗ | 1.0 ✗ (now too *few*) |
| SO coupling | 50–70% | 35–86% ✗ | 100% ✗ (too deterministic) |
| | | **2/5** | **6/8** |

**The events now qualify as sleep spindles**: 0.5 s trains at 13.7 Hz built from
genuine rebound bursts. The two remaining failures have flipped direction and
are tuning rather than mechanism — spindles are now slightly too **rare**
(1/min vs 2–8) and too **perfectly** locked to the SO (100% vs 50–70%, because
every spindle is initiated by the deterministic trigger). Both should yield to
raising thalamic excitability slightly and adding jitter/spontaneous initiation.

**Cost.** AdEx routes E/I by weight **sign**, so it has **no separate GABA_B
receptor**; the slow GABA_B component of RE→TC (Mushtaq Table 3) is not
represented in the AdEx configuration, whereas `ht_neuron` does model it. Both
models remain selectable (`config/network_auditory_adex.yaml`,
`config/network_auditory_hh.yaml`).

## 5.7 Final result: 10/10 criteria

Adding **spontaneous, non-SO-locked spindle initiation** closes the last two
gaps. The review (sect. VI.E) states only 50–70% of spindles are SO-coupled —
i.e. **30–50% occur with no SO co-detected** — so a purely deterministic,
SO-locked trigger is itself unphysiological. A low-rate Poisson drive onto nRT
(`SleepParams.spontaneous_re_rate`) lets the reticular nucleus initiate spindles
on its own, raising density *and* breaking the deterministic SO locking.

**Calibration** (AdEx, 40 TC / 40 RE, 150–200 s). The usable window is narrow:

| spontaneous rate | density | duration | note |
|---|---|---|---|
| 0 Hz | 0.8 /min | 0.59 s | too few |
| **4 Hz** | **2.0 /min** | **0.68 s** | chosen |
| 6 Hz | 1.5 /min | 0.61 s | declining |
| 8 Hz | 0.6 /min | 0.77 s | declining |
| 12 Hz | 0.0 /min | — | spindles collapse (nRT driven out of burst mode) |

**Validation, AdEx + 4 Hz spontaneous, 300 s:**

| criterion | paper | measured | |
|---|---|---|---|
| intra-spindle frequency | 10–15 Hz | 13.7 Hz | ✓ |
| duration | 0.5–3 s | 0.66 s | ✓ |
| density | 2–8 /min | 2.4 /min | ✓ |
| **inter-spindle interval** | **5–10 s** | **6.6 s** | **✓** |
| SO coupling | 50–70% | 50% | ✓ |
| infraslow clustering | ~0.02 Hz | 0.020 Hz (ratio 432) | ✓ |
| RE spikes per burst | 2 to >10 (>100 Hz) | 5.24 | ✓ |
| TC spikes per burst | 2–6 (>100 Hz) | 4.27 | ✓ |
| RE V_m < −55 mV | required | 97% | ✓ |
| TC V_m < −65 mV | required | 36% | ✓ |
| | | **10/10** | |

Notably the **5–10 s refractory period** — which resisted three direct attempts
(deeper I_h, trigger spacing, gated mode) — falls out for free once the cells
genuinely burst: after a rebound burst the AdEx adaptation current `w` keeps the
relay refractory for seconds. It was never a drive-timing problem.

**Statistical caveat.** Spindle density is low by construction (2–8/min), so
short runs contain very few events and the SO-coupling percentage is unstable:
across 150–200 s runs it ranged 0–80% with n = 2–5. **Validate over ≥300 s**
(n ≈ 10+); the earlier "100% SO coupling" was computed from roughly one event
and was meaningless.

## 5.8 Correction: ht_neuron CAN rebound-burst (bio-plausibility restored)

**Correcting sect. 5.5.1.** That section concluded ht_neuron's I_T "does not
produce a supra-threshold rebound". **That was wrong**, and the error propagated
into the decision to adopt AdEx. Re-tested properly, ht_neuron produces textbook
rebound bursts across a wide regime:

| g_KL | g_peak_T | rebound |
|---|---|---|
| 0.5 | 3 | 6 spikes @ 107 Hz |
| 0.6 | 6 | 5 spikes @ 167 Hz |
| 1.0 | 3 | 2 spikes @ 118 Hz |
| 1.0 | 10 | 4 spikes @ 210 Hz |

The original calibration failed because it used too-shallow hyperpolarisation
(-60 vs -100 pA), too short a pulse (I_T de-inactivation is slow), and
g_peak_T = 1-3. Full +30 mV action potentials, 2-6 per burst, exactly as the
review describes.

**Why it still failed synaptically — a real structural bind.** Reaching the I_T
de-inactivation window needs V < ~-80 mV, but:

- **GABA_A** is fast (7 ms decay, quick release -> good for rebound) but its
  reversal is a hard floor: with E_GABA_A = -70 mV it cannot take the cell below
  -70 mV, and at high g_KL the cell rests *below* -70 so GABA_A becomes
  shunting/depolarising. Measured: Vmin -69.9 mV, no rebound.
- **GABA_B** reaches -89 mV, but its 200 ms decay then **clamps** the membrane
  at -72...-77 mV long after the burst, so I_T re-inactivates before the cell
  can rebound. Measured: Vmin -89.5 mV, still zero rebound spikes.

**Resolution (paper-supported).** The review (sect. V.B.3) states TC cells "show
a hyperpolarized reversal potential for Cl-, which leads to large inhibitory
postsynaptic potentials". Setting `HHParams.E_rev_GABA_A_TC = -85 mV` gives deep
inhibition through the FAST GABA_A conductance, so the cell is released quickly
and I_T fires: **2-3 rebound spikes at 87-126 Hz** from a synaptic RE volley.

**Network result** (ht_neuron, 40 TC / 40 RE, 300 s): **5/9**, up from 2/5.
Duration 0.53 s PASS, intra-spindle frequency 11.3 Hz PASS, infraslow 0.020 Hz
PASS, V_m ranges PASS. Still failing: in-network bursting (RE 1.45, TC 1.11
spikes/burst) and density 0.2/min -- the single-cell fix does not fully
transfer, because in-network RE volleys are not synchronised enough to deliver
the deep, fast IPSP the single-cell test used.

### Model choice: a genuine trade-off

| | `ht_neuron` (HH) | `aeif_cond_exp` (AdEx) |
|---|---|---|
| criteria passed | 5/9 | **10/10** |
| biophysics | **conductance-based: I_T, I_h, I_NaP, I_KNa** | phenomenological adaptation `w` |
| rebound mechanism | **real T-type Ca current** | fitted 2-variable dynamics |
| GABA_B | **yes** | no (weight-sign routing) |
| in-network bursting | not yet | yes |

AdEx passes more criteria but is an integrate-and-fire model whose "rebound
burst" is a curve fitted to the phenomenon rather than an ion current; ht_neuron
is substantially more bio-plausible and is the right target for mechanistic
work. **Remaining gap for HH: in-network RE burst synchrony** (the review's
sect. V.C mechanisms -- gap junctions, topographic slabs -- which NEST/ht_neuron
partly cannot express).

## 6. References

- Fernandez LMJ, Lüthi A. *Sleep Spindles: Mechanisms and Functions.*
  Physiol Rev 100: 805–868, 2020.
  [doi:10.1152/physrev.00042.2018](https://doi.org/10.1152/physrev.00042.2018)
- Mushtaq M, Marshall L, ul Haq R, Martinez T. PLOS ONE 19(6):e0306218, 2024 —
  source of the current synaptic parameters (see the README).
- Hill S, Tononi G. *J Neurophysiol* 93: 1671–1698, 2005 — the `ht_neuron` model.
- See also [`model_equations.md`](model_equations.md) for the implemented
  neuron and synapse equations.
