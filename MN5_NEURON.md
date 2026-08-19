# MN5 — NEURON MPI run

The NEURON path (HH + SK2 + Ca-dependent I_h). This supersedes the NEST
manifest in [`MN5_UPLOAD.md`](MN5_UPLOAD.md) for scale-out: NEST's point
neurons have no Ca²⁺ pool and no Ca²⁺-dependent SK2, so they cannot express the
mechanism under test. Keep the NEST files only if you also want the reference
results.

## Files to upload

Only 10 files — `tc_mpi.py` imports `tc_neuron` and nothing else from the repo.

| path | why |
|------|-----|
| `run_nrn.sh` | SLURM script (100 tasks) |
| `neuron/tc_mpi.py` | the MPI network + benchmark |
| `neuron/tc_neuron.py` | `TCCell` / `RECell` — the cell definitions |
| `neuron/mod/itd.mod` | Destexhe T-current (the burst generator) |
| `neuron/mod/hh2.mod` | Traub–Miles Na⁺/K⁺ spikes |
| `neuron/mod/cad.mod` | submembrane Ca²⁺ pool |
| `neuron/mod/sk2.mod` | SK2 Ca²⁺-activated K⁺ |
| `neuron/mod/ihca.mod` | Ca²⁺-dependent I_h |
| `neuron/mod/gapmpi.mod` | **cross-rank** gap junction (RANGE `vgap`) |
| `neuron/mod/gap.mod` | serial gap junction (`nrnivmodl` compiles the whole dir) |

`requirements.txt` is optional — the run needs only `neuron` + `numpy`, plus
`scipy` for the Welch spectral estimator (it falls back to a Hann-windowed FFT
if scipy is absent).

**Not needed:** any `tc_network*.py`, `config/*.yaml`, `tc_validate.py`,
`tc_run.py`, or the NEST scripts. `tc_mpi.py` takes its parameters as CLI flags.

## scp

```bash
ssh USER@glogin1.bsc.es 'mkdir -p ~/audCC/neuron/mod'
scp run_nrn.sh USER@glogin1.bsc.es:~/audCC/
scp neuron/tc_mpi.py neuron/tc_neuron.py USER@glogin1.bsc.es:~/audCC/neuron/
scp neuron/mod/*.mod USER@glogin1.bsc.es:~/audCC/neuron/mod/
```

## On MN5

**Do not pip install anything** — MN5 does not permit it, and NEURON is already
provided as a site module. `run_nrn.sh` loads modules rather than building a
venv.

Module order matters (Lmod): `python/3.12.1` has a hard dependency on `intel`,
so `module load python` on its own fails with *"Cannot load module
python/3.12.1 without these module(s) loaded: intel"*. The script loads `intel`
first.

Find the NEURON module name once:

```bash
module spider neuron
```

If it is not called `neuron`, pass the real name on the sbatch line
(`NEURON_MODULE=...`); `INTEL_MODULE` and `PYTHON_MODULE` override the other two
the same way.

Compile the mechanisms once on the login node (keeps a serial build out of your
100-CPU allocation and surfaces compiler errors immediately):

```bash
module load intel python neuron
cd ~/audCC/neuron && nrnivmodl mod
```

`scipy` is **optional** — the spectral estimator falls back to a Hann-windowed
FFT, verified to agree with the scipy path to **0.01 Hz**. Only `neuron` and
`numpy` are required.

Edit `run_nrn.sh` to set `--account` / `--qos` for your project (both are
commented out at the top), then:

### 1. Benchmark first

```bash
sbatch --export=ALL,BENCH=1 run_nrn.sh
```

Times 100 / 500 / 2000 / 5000 cells at 1 s each. This replaces the laptop
`N^1.67` extrapolation — which was measured thread-parallel on one machine and
is **not** a measurement of NEURON-MPI — and tells you whether 5k HH cells is
affordable, and whether 100 ranks is the right request, before you spend a
production run finding out.

### 2. Production, once the benchmark justifies the size

```bash
sbatch --export=ALL,NTC=2500,NRE=2500,TSTOP=200000 run_nrn.sh
```

`TSTOP=200000` (200 s) is what the spindle density / ISI statistics need.
Other knobs: `GH` (default `4e-4`), `GSK` (`5e-5`), `TAG`.

Results land in `out/<TAG>.npz` (spike times + gids + metadata):

```bash
scp USER@glogin1.bsc.es:'~/audCC/out/*.npz' out/
```

## Before trusting a production number

Two open items, both recorded in [`neuron/README.md`](neuron/README.md):

- **`GH=4e-4` is ~20× Destexhe's published `ghbar = 2×10⁻⁵`.** The 13 Hz spindle
  rhythm depends on it, and that departure needs justification before it counts
  as a result.
- **The `gh_tc` response is non-monotonic** — `2e-4` gives a broken, irregular
  3 Hz regime between the 9 Hz and 13 Hz points. Unexplained.

The mechanism and topology themselves are verified: the network is
**rank-count independent** (13.69 Hz, IQR 35 ms, n=114, bit-identical at 1, 2
and 4 ranks) and consistent with the serial model (13.41 Hz, IQR 31 ms).

## Corticothalamic scale-out (Aug 2026): `ctx_thalamus_mpi.py`

The thalamus-only path above tests the RE↔TC loop in isolation. Once cortex
is real (`neuron/cortex_neuron.py`, `neuron/ctx_thalamus_network.py` — L5
intrinsic-bursting cells generating the column's own slow oscillation, closed
onto TC/RE via L6), the loop needs to run at a bio-plausible cortical column
size, not the 10+10 (thalamus) + ~40-cell (cortex) laptop demo. That is what
`neuron/ctx_thalamus_mpi.py` is for: the same `ParallelContext` / fixed-
convergence pattern as `tc_mpi.py`, extended with the full L4/L2-3/L5/L6
column.

### Files to upload (adds 4 to the thalamus-only list)

| path | why |
|------|-----|
| `run_ctx_nrn.sh` | SLURM script (companion to `run_nrn.sh`) |
| `neuron/ctx_thalamus_mpi.py` | the MPI corticothalamic network + benchmark |
| `neuron/cortex_neuron.py` | `PYCell` / `PYCellIB` / `FSCell` — cortical cell definitions |
| `neuron/mod/ical.mod` | HVA Ca²⁺ (cortical SK2 Ca²⁺ source) |
| `neuron/mod/inap.mod` | persistent Na⁺ (L5 `PYCellIB` intrinsic bursting) |

Everything else (`tc_neuron.py`, `it`/`it2`/`hh2`/`cad`/`sk2`/`ihca`/`gap`/
`gapmpi` mod files) is already on the thalamus-only list above — `nrnivmodl
mod` compiles the whole `mod/` directory in one pass regardless.

### Population sizes

Default sizes sum to **3050** cells, the DECLARED size in
`config/network_auditory_mn5_5k.yaml` (the NEST reference model reaches
~5010 after its own RS/FRB/TuftRS-TuftIB/Basket-LTS-Axoaxonic subtype
splits, which this NEURON port does not reproduce 1:1):

```
thalamus  TC 210 / RE 55
L4        E 640 / I 160
L2/3      E 640 / I 160
L5        E 530 / I 130   (half PYCellIB, half PYCell)
L6        E 420 / I 105
```

`--scale 1.65` multiplies every population by 1.65×, landing at ~5030 cells
— closer to the NEST model's actual (post-split) size, if you want the
literal "5k" figure rather than the declared 3050.

### 1. Benchmark first

```bash
sbatch --export=ALL,BENCH=1 run_ctx_nrn.sh
```

Times 0.1×/0.5×/1×/1.65× the default sizes (≈300/1500/3050/5030 cells), 1 s
each — same reasoning as `run_nrn.sh --bench`: measure before committing a
200 s production run.

### 2. Production

```bash
sbatch --export=ALL,SCALE=1.0,TSTOP=200000 run_ctx_nrn.sh
```

`SCALE=1.0` → the declared 3050-cell reference size; `SCALE=1.65` → ~5030.
Other knobs: `CONV` (fixed convergence, default 100 — see `tc_mpi.py`'s note
on why fixed convergence and not the serial model's `frac`-based fan-in),
`TAG`.

### ⚠️ Status: locally validated at 1 rank only — NOT yet run on MN5, and
### multi-rank correctness is UNTESTED here

This machine has no MPI runtime NEURON can link against (`libmpi.dylib` not
found — `mpirun -n 4 nrniv -python -mpi ...` fails at `nrnmpi_init`), so only
`nhost=1` could be exercised locally. What IS verified at `nhost=1`:

- The network builds without error at the reference gid layout and produces
  sensible dynamics — L5's `PYCellIB` population self-generates its slow
  oscillation and the signal correctly propagates L5→L6→(TC, RE): a
  61-cell smoke run (`--scale 0.02 --tstop 3000`) shows TC/RE/L6 all spiking,
  not just L5.
- `--bench` runs and reports a sane `s per 1k cells` figure at small scale.

What is **NOT** verified, because it needs ≥2 real MPI ranks:

- Cross-rank `gid_connect` (every inter-population projection).
- Cross-rank gap junctions (`pc.setup_transfer` / `GapMPI`) between RE cells.
- Rank-count independence of the result (the thalamus-only path's own
  bit-identical-across-ranks claim, above, took a real multi-rank run to
  establish — this module inherits the same code pattern but has not had the
  same check).

Treat this as ready to **attempt** on MN5, not as validated. Run the
benchmark at a small `--bench-scales` (e.g. `0.02,0.05`) with `--ntasks 2` and
`--ntasks 4` FIRST and diff the spike output before trusting a 200 s / 100-rank
production run — the same rank-count-independence check the thalamus-only
path already passed.

### The `gh_tc = 4e-4` caveat above does NOT apply to this network

`ctx_thalamus_mpi.py`'s `tc`/`re` cells hard-code `gh=0.0` and
`g_tc_re=0.011` (matching `ctx_thalamus_network.py`'s current defaults, see
`neuron/README.md` "corticothalamic loop" section), not the thalamus-only
model's `gh_tc=4e-4` / `g_tc_re=0.004`. That thalamus-only setting was found
to make an isolated `TCCell` fire spontaneously at ~103 Hz with zero synaptic
input; it is not carried over here. The combined model's own open item is
different: RE↔TC produces 1–5 oscillatory cycles at ~7 Hz per SO event, not
yet a consistent 10–15 Hz spindle train — see `neuron/README.md` for the full
account. **Whether a bio-plausible cell count changes this is the open
question an MN5 run is meant to answer** — small-N quantization (only 10 TC +
10 RE cells drawing from a handful of possible wiring motifs) is one
candidate explanation, tested at small scale in
`neuron/ctx_thalamus_network.py`'s new `het=`/`delay_jitter=` heterogeneity
parameters before committing to this larger run (see `neuron/README.md`
"Heterogeneity" section for that result).
