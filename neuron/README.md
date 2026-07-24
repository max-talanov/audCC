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
| `mod/cav3.mod` | Ca_v3.1 / Ca_v3.3 low-threshold T-current (Huguenard & McCormick kinetics) |
| `mod/hh2.mod` | Traub-Miles fast Na⁺ / K⁺; fires repetitively on a plateau (no depol. block) |
| `mod/cad.mod` | submembrane Ca²⁺ pool (feeds SK2) — experimental |
| `mod/sk2.mod` | SK2 Ca²⁺-activated K⁺ burst terminator — experimental |
| `tc_neuron.py` | two-compartment TC relay cell + single-cell rebound-burst demo |
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
- **Multi-spike rebound BURST achieved.** Swapping the built-in `hh` soma for
  Traub–Miles `hh2` kinetics **solves the depolarisation block**: the soma now
  fires a genuine **repetitive rebound burst** on the T-current plateau (21 / 26
  / 31 spikes at gCaT = 0.01 / 0.02 / 0.05, and **0** at gCaT = 0 — the burst is
  T-current-driven, not an artefact). This is the mechanism no NEST point model
  (`ht_neuron`, AdEx) could produce.

- **SK2 burst terminator is now functional** (`mod/sk2.mod` + `mod/cad.mod`,
  opt-in via `gsk > 0`). The earlier cad↔cav3 feedback bug is **fixed**: `cad`
  now accumulates into a **private ion species** (`sk`, valence 2) that only SK2
  reads, so the submembrane pool no longer perturbs the T-current's Ca²⁺
  reversal. With this, the rebound burst survives (`gsk = 0` → ~13 spikes) and
  SK2 **controls its length** — raising `gsk` shortens it (→ 2 spikes). `depth`
  (Ca accumulation rate) and `kd` set when SK crosses threshold.

**Next**
- **A *tight* physiological 2–6 spike burst.** SK2 shortens the burst but not yet
  into a clean 2–6 spike / ~30 ms packet: it currently yields a couple of spikes
  spread across the LTS. Root cause is upstream of SK2 — the low-threshold
  Ca²⁺ **plateau is too long (~250 ms)** because somatic spikes back-propagate
  and *re-prime* the dendritic T-current. Fixing this needs a richer relay-cell
  model (proper dendritic geometry + additional currents, as in Destexhe's full
  TC cell), not more SK2 tuning. This is the honest remaining gap.
- Reticular (RE) cell with Ca_v3.3 + SK2 and **gap junctions** (the synchrony
  mechanism NEST could not provide).
- Port the network topology; validate against the same 10-criteria harness.
- The NEST model stays the working reference throughout; consider keeping it for
  MareNostrum 5 scale-out (NEST is the better large-network tool).
