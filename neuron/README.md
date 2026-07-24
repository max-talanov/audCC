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
| `tc_neuron.py` | two-compartment TC relay cell + single-cell rebound demo |
| `arm64/` (git-ignored) | compiled mechanism, from `nrnivmodl mod` |

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
  hyperpolarisation** (h → 0.99) and drives a strong **inward Ca current that
  scales with gCaT** (−0.7 to −1.3 mA/cm²) on release — the exact mechanism the
  NEST point models could not express. It produces a post-inhibitory rebound to
  spike threshold.

**Next**
- **Multi-spike rebound burst.** The T-current fires, but the HH soma emits only
  one Na⁺ spike and then enters depolarisation block on the Ca²⁺ plateau. A full
  2–6 spike burst needs proper delayed-rectifier / A-type K⁺ kinetics in the
  soma (HH's Na/K do not support bursting here) — a dedicated soma `.mod`.
- Reticular (RE) cell with Ca_v3.3 + SK2 and **gap junctions** (the synchrony
  mechanism NEST could not provide).
- Port the network topology; validate against the same 10-criteria harness.
- The NEST model stays the working reference throughout; consider keeping it for
  MareNostrum 5 scale-out (NEST is the better large-network tool).
