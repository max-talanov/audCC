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
