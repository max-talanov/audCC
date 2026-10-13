# Plan: a richer auditory thalamocortical architecture

Status (2026-10-12): **plan only, nothing implemented.** Compares the
current model with the Traub et al. (2005) single-column thalamocortical
model (14 cell types; the picture in `max-talanov/cc`, `CC.png`) and with
known auditory-cortex cell types, and stages the additions. Runs alongside
`PLAN-spindels.md` and `PLAN-auditory-input.md`, and takes over their
architecture items: TC → L4I (auditory D3) and depressing / facilitating
synapses (auditory B2).

## Goal

Add the cell types and connection classes that matter for this project's
questions:
- the sleep slow oscillation (SO);
- spindles;
- how sounds are relayed and gated in sleep and wake;
- later, gamma / 40 Hz following and stimulus-specific adaptation (SSA).

Every addition sits behind a flag with the old behaviour as default, and
must keep the full-scale sleep model's results (`PLAN-spindels.md` Stage 3,
`PLAN-auditory-input.md` C1/C2) or explain why they change.

## The current architecture

Full scale (`SCALE=1.65`, `SPINDLE=stage2 --cx-noise-rate 200`, K_Na), counted
by building the network; diagrams in `figures/architecture_cortex.png`,
`figures/architecture_thalamocortical.png` (`neuron/plot_architecture.py`).

- **Neurons: 5,031** in 10 populations. There are 4 cell models: `PYCell`
  (RS), `PYCellIB`, `FSCell`, plus `TCCell` / `RECell`.
  - L2/3 E 1056 / I 264;
  - L4 E 1056 / I 264;
  - L5 E 874 (437 IB + 437 RS) / I 214;
  - L6 E 693 / I 173;
  - TC 346, RE 91.
- **Chemical synapses: 1,158,676.** Fixed convergence of 100 random sources
  per target and projection; the thalamus is wired locally on a ring.

  | Group | Synapses |
  |---|---|
  | feed-forward E → E: L4 → L2/3, L2/3 → L5, L4 → L5, L5 → L6, L5 → L2/3 | 455,300 |
  | L5 recurrent NMDA | 87,400 |
  | E → I within layer | 91,500 |
  | I → E within layer | 367,900 |
  | TC → L4E | 105,600 |
  | L6 → TC / RE | 43,700 |
  | TC ↔ RE, RE ↔ RE | 7,276 |

- **Gap junctions:** RE (715 halves), L5 IB (3,481).
- **External input:** background Poisson input on cortical E cells; IC → TC
  tones.
- **Interneurons are 19.9% of cortical neurons**, which is realistic, but
  all of them are one fast-spiking basket-like type.

## Gap analysis

### Against Traub et al. 2005

Diagram: `figures/architecture_traub.png` (`neuron/plot_architecture.py`). Counts are from the
paper's Fig. 1 legend: 14 types, 3,560 neurons. Of the 14, 9 are present (one as a pyramidal
stand-in) and 5 are missing: SupFRB, SupAxAx, SupLTS, DeepAxAx, DeepLTS.

| Traub type | Ours | Status |
|---|---|---|
| SyppyrRS (L2/3 RS pyramid) | L2/3 E | ✅, but no L2/3 → L2/3 recurrence |
| SyppyrFRB (L2/3 fast-rhythmic-bursting) | — | ❌ |
| Bask23 (L2/3 basket) | L2/3 I | ✅ |
| Axax23 (L2/3 axo-axonic) | — | ❌ |
| LTS23 (L2/3 low-threshold-spiking) | — | ❌ |
| Spinstel4 (L4 spiny stellate) | L4 E (`PYCell`) | ⚠️ pyramidal stand-in; its only input is TC |
| LTS4 | — | in the cc model's `CC.png`, but not in Traub 2005: his layer 4 has no interneurons (its inhibition comes from the deep ones). Stage 2 adds it as an auditory-cortex SST population |
| TuftRS5 / TuftIB5 | L5 RS / IB | ✅ |
| Bask56 | L5 I, L6 I | ✅ |
| Axax56 | — | ❌ |
| LTS56 | — | ❌ |
| NontuftRS6 | L6 E | ✅ |
| TCR / nRT | TC / RE | ✅ |

Missing connection classes:
1. E → E recurrence within L2/3, L4, L6 (only L5 has it);
2. any I → I synapse;
3. cross-layer inhibition (e.g. LTS23 / Bask23 → L5 tufts);
4. TC → interneurons (feed-forward inhibition) and TC → L5 / L6;
5. dendritic gap junctions within interneuron types, and axonal gap
   junctions between pyramids;
6. GABA_B.

### Auditory-cortex types beyond Traub

| Type | Relevance here |
|---|---|
| SST / Martinotti (Traub's LTS) | Dendritic, facilitating inhibition. PV and SST cells shape SSA in auditory cortex differently (Natan et al. 2015). SST cells promote sleep slow waves (Funk et al. 2017) |
| VIP (disinhibition: VIP → SST) | State and arousal switch; mainly for wake (auditory Stage B) |
| L1 neurogliaform cells | Slow GABA_A / GABA_B, apical-tuft gating |
| L5 IT vs PT; PT corticocollicular | PT cells project back to the inferior colliculus: a feedback loop onto our IC input |
| L6 corticothalamic vs L6 IT | Only corticothalamic cells project to TC / RE; now all L6 E do |
| Non-lemniscal thalamus (MGBd / MGBm) → L1 / L5 | Diffuse thalamic drive (out of scope in the spindle plan) |

## Design rules

- **Gids stay stable.** Populations are contiguous gid blocks in
  `DEFAULT_SIZES` order, and every per-cell parameter and wiring draw is a
  function of gid. New populations are **appended after `l6i` with default
  size 0**, so with all flags off every existing gid, cell and synapse is
  unchanged, and old seeds reproduce. A gid-identical regression check (as
  in auditory Stage A, `aud_analyze.py --compare`) is the gate for each
  stage.
- **Interneuron classes split the existing interneuron budget** rather than
  adding cells, keeping ~20% inhibitory. Target split at scale 1.65, a
  starting point from general mouse-cortex proportions (Tremblay et al.
  2016; PV ~40%, SST ~30%, VIP / 5HT3aR ~30%, VIP mostly superficial). Check
  against auditory-cortex counts before freezing:

  | Layer | I now | PV (FS) | SST (LTS) | VIP |
  |---|---|---|---|---|
  | L2/3 | 264 | 106 | 79 | 79 |
  | L4 | 264 | 172 | 66 | 26 |
  | L5 | 214 | 107 | 86 | 21 |
  | L6 | 173 | 87 | 69 | 17 |

  A flag (`--arch traub` / per-class fractions) moves cells from the FS
  population into the new classes. The FS gid blocks keep their positions,
  and the moved cells live in the appended blocks.
- **Analysis must follow** (`ctx_analyze.py`, `stage2_report.py`,
  `aud_phase.py`, the LFP proxy): new populations get counted and plotted.
  The UP-state detector keeps using L5E + L6E.

## Stage 1 — new connections, no new cell types (local, then MN5)

- [ ] **TC → L4 FS** (feed-forward inhibition; auditory D3) and **TC → L5 /
      L6 E**: `--g-tc-l4i`, `--g-tc-l5`, `--g-tc-l6`.
- [ ] **E → E recurrence in L2/3, L4, L6**, AMPA plus optional NMDA as in L5:
      `--g-l23-rec`, `--g-l4-rec`, `--g-l6-rec`.
- [ ] **I → I**: FS → FS within each layer (`--g-i-i`).
- [ ] Regression with the flags off: gid-identical.
- [ ] With the flags on, re-check at scale 0.1 (5 seeds × 20 s) and on MN5
      (200 s):
  - SO 0.5–1 Hz with CV ≥ 0.3;
  - spindles ≥ 6 cycles, 10–15 Hz, ~0.5 s, nested in UP states;
  - and with tones, C1/C2.

**Expected / to test:** TC → FS should sharpen tone responses and may reduce
the all-or-none evoked spindles (`PLAN-auditory-input.md`, "C1/C2 final").
Recurrence in L2/3 and L4 lets UP states involve the upper layers, which are
now driven only from L4 / L5.

**Done when:** each new connection's effect on the SO, spindles and C1/C2 is
measured, the useful ones are adopted into the sleep preset, and the flags
off reproduce the old runs.

## Stage 2 — SST / LTS interneurons (local, then MN5)

- [ ] Cell model `LTSCell`: hh2 + T-type Ca²⁺ (low-threshold rebound) + spike
      adaptation (SK2 / M-current), calibrated to LTS firing (adapting
      trains, rebound bursts).
- [ ] Populations `l23s`, `l4s`, `l56s` (appended gids), sized from the
      table above.
- [ ] Synapses:
  - SST → E on the **dendrite** (`synsec`), slower GABA_A, **facilitating**.
    This needs a short-term-plasticity synapse, `mod/tmgsyn.mod`
    (Tsodyks–Markram), shared with auditory B2;
  - E → SST facilitating;
  - SST → FS;
  - **cross-layer** L2/3 SST → L5 E tufts (Traub's LTS23 → TuftIB / RS).
- [ ] Tests:
  - DOWN-state onset and SO regularity (SST and slow waves);
  - L5 IB bursting;
  - tone responses in UP vs DOWN;
  - later, SSA with an oddball protocol (auditory B2).

**Done when:** SST populations fire as LTS cells, the sleep model still
meets the Stage 3 criteria, and their effect on the SO and on tone
responses is documented.

## Stage 3 — FRB (chattering) pyramids in L2/3

- [ ] Cell model: `PYCell` + persistent Na⁺ / fast AHP tuned to fast rhythmic
      bursting (bursts at 30–50 Hz on depolarisation).
- [ ] Population `l23frb` (part of L2/3 E, appended), with FRB → FRB and FRB ↔
      RS recurrence.
- [ ] Test, in wake (needs auditory Stage B): gamma and the 40 Hz click-train
      response (auditory B2).

## Stage 4 — VIP interneurons (with auditory Stage B)

- [ ] Cell model: an adapting, irregular-spiking interneuron.
- [ ] VIP → SST (disinhibition), weak VIP → FS; driven by the wake state
      (neuromodulation in `brain_state.py`).
- [ ] Test: wake / NREM contrast of upper-layer tone responses.

## Stage 5 — output-type splits and auditory feedback

- [ ] L5 E → L5 IT / L5 PT; L6 E → L6 corticothalamic / L6 IT. Only
      corticothalamic cells project to TC / RE; PT projects out.
- [ ] **Corticocollicular feedback**: L5 PT → IC fibres (a gain or
      rate-modulating input on `auditory_input.py`).

## Out of scope for now

- Axo-axonic cells (Traub's Axax23 / 56) and axonal gap junctions between
  pyramids. They mainly matter for Traub's very fast oscillations and
  epileptiform bursts.
- L1 neurogliaform cells, GABA_B (revisit if the SO DOWN states need slower
  inhibition).
- Non-lemniscal thalamus (MGBd / MGBm).
- Detailed multi-compartment morphologies (tufts as separate compartments).
  The cells stay at 1–2 compartments, with "dendritic" synapses on `synsec`.

## Dependencies

| This plan | Needs | Gives |
|---|---|---|
| Stage 1 | — | auditory D3 (TC → L4I) |
| Stage 2 | `tmgsyn.mod` (shared with auditory B2) | auditory B2 (SSA), spindle / SO checks |
| Stage 3, 4 | auditory Stage B (a calibrated wake state) | auditory B2 (40 Hz), Stage B contrast |
| Stage 5 | auditory D1 (IC input exists) | closed auditory loop |

## Risks

- **The sleep model is finely tuned.** New inhibition or recurrence can
  break the SO (0.77 Hz, CV 0.38) or the spindles. Each stage re-runs the
  Stage 3 checks, and parameters are re-balanced one at a time, as in the
  spindle plan.
- **Proportions and connection probabilities from mouse cortex, not
  auditory cortex specifically.** They are treated as starting points and
  documented as such.
- **More synapses, more cost.** Stage 1 adds ~3 × 100 × (L2/3 + L4 + L6 E) ≈
  0.3–0.9 M synapses. That is fine on MN5 (a 5k-cell run is ~10 min / 200 s),
  but slower for local scale-0.1 sweeps.

## References

- Traub RD et al. (2005) Single-column thalamocortical network model
  exhibiting gamma oscillations, sleep spindles, and epileptogenic bursts.
  J Neurophysiol 93: 2194–2232.
- Tremblay R, Lee S, Rudy B (2016) GABAergic interneurons in the neocortex:
  from cellular properties to circuits. Neuron 91: 260–292.
- Natan RG et al. (2015) Complementary control of sensory adaptation by two
  types of cortical interneurons. eLife 4: e09868.
- Funk CM et al. (2017) Role of somatostatin-positive cortical interneurons
  in the generation of sleep slow waves. J Neurosci 37: 9132–9148.
- Pi H-J et al. (2013) Cortical interneurons that specialize in
  disinhibitory control. Nature 503: 521–524.
- Fernandez LMJ, Lüthi A (2020) Sleep spindles: mechanisms and functions.
  Physiol Rev 100: 805–868.
