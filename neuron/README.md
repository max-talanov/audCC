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
| `mod/cad.mod` | submembrane Ca²⁺ pool in a private `sk` ion (feeds SK2) |
| `mod/sk2.mod` | SK2 Ca²⁺-activated K⁺ burst terminator |
| `mod/gap.mod` | electrical (gap-junction) coupling between TRN cells (connexin-36) |
| `tc_neuron.py` | `TCCell` (relay), `RECell` (reticular), `gap_junction()` + demos |
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

**Next**
- Port the network topology (TC↔RE loop, corticothalamic drive); validate
  against the same 10-criteria harness (`tc_validate.py`, simulator-agnostic).
- The NEST model stays the working reference throughout; consider keeping it for
  MareNostrum 5 scale-out (NEST is the better large-network tool).
