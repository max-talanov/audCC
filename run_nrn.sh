#!/bin/bash
#SBATCH --job-name=audcc_nrn
#SBATCH --output=slurmout_%j.txt
#SBATCH --error=slurmerr_%j.txt
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --qos=gp_debug
#
# MareNostrum 5 -- MPI-parallel NEURON thalamic network (HH + SK2 + Ca-dependent
# I_h). This replaces the NEST path: NEST's point neurons cannot express the
# SK2 / low-threshold-Ca balance, so they cannot test the hypothesis.
#
#   sbatch --export=ALL,BENCH=1 run_nrn.sh                 # scaling benchmark
#   sbatch --export=ALL,NTC=2500,NRE=2500,TSTOP=200000 run_nrn.sh   # production
#
# RUN THE BENCHMARK FIRST. It replaces a laptop N^1.67 extrapolation with a real
# measurement and tells you whether 5k cells is affordable before you spend a
# 200 s production run on it.

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$WORKDIR"
[ -f neuron/tc_mpi.py ] || { echo "ERROR: run sbatch from the repo root (no neuron/tc_mpi.py in $WORKDIR)"; exit 1; }

# --- modules ---------------------------------------------------------------
# ORDER MATTERS on MN5 (Lmod): python/3.12.1 declares a hard dependency on
# intel, so `module load python` alone fails with
#   "Cannot load module python/3.12.1 without these module(s) loaded: intel"
# Nothing is pip-installed here -- MN5 does not permit it. NEURON comes from a
# site module. Override any of these from the sbatch --export line if the names
# differ; find them with `module avail neuron` / `module spider neuron`.
# --- environment ------------------------------------------------------------
# This script loads NO modules, matching run.sh. The environment is set up
# OUTSIDE it (e.g. `module load miniforge`, which provides python + NEURON).
#
# Earlier versions loaded modules here and failed three times over, each time
# for a different reason: `module purge` stripped MN5's default mkl/impi, then
# python/3.12.1 wanted intel, then mkl -- and finally python/3.12.1 CONFLICTS
# with an already-loaded miniforge ("Cannot load module python/3.12.1 because
# these module(s) are loaded: miniforge"). Whatever module set you use, load it
# before sbatch and leave this script alone.
#
# If you ever need a module here, uncomment and adjust:
# module load miniforge
# module load neuron

PY="${PY:-}"
if [ -z "$PY" ]; then
    if [ -x ./.venv-neuron/bin/python ]; then PY=./.venv-neuron/bin/python
    else PY="$(command -v python3 || command -v python)"; fi
fi
[ -n "$PY" ] && [ -x "$PY" ] || { echo "ERROR: no python found (set PY=... )"; exit 1; }
echo "python : $PY"

"$PY" -c "import neuron" 2>/dev/null || {
    echo "ERROR: '$PY' cannot import neuron."
    echo "       NEURON is provided by a module on MN5 -- do NOT pip install."
    echo "       Try:  module spider neuron   then resubmit with NEURON_MODULE=<name>"
    exit 1
}
"$PY" -c "import numpy" 2>/dev/null || { echo "ERROR: numpy unavailable"; exit 1; }
# scipy is OPTIONAL -- the spectral estimator falls back to a Hann-windowed FFT.
"$PY" -c "import scipy" 2>/dev/null || echo "note: no scipy; using the FFT fallback"

# --- compile the NMODL mechanisms (itd, hh2, cad, sk2, ihca, gap, gapmpi) ----
cd neuron
if [ ! -f x86_64/libnrnmech.so ]; then
    echo "== compiling NMODL mechanisms =="
    command -v nrnivmodl >/dev/null || { echo "ERROR: nrnivmodl not on PATH -- is the NEURON module loaded?"; exit 1; }
    nrnivmodl mod
fi
cd "$WORKDIR"

NTC="${NTC:-2500}"
NRE="${NRE:-2500}"
TSTOP="${TSTOP:-200000}"
GH="${GH:-4e-4}"
GSK="${GSK:-5e-5}"
TAG="${TAG:-nrn_${SLURM_JOB_ID:-local}}"
mkdir -p out

echo "== NEURON MPI: ${SLURM_NTASKS:-1} ranks =="

if [ "${BENCH:-0}" = "1" ]; then
    # Scaling benchmark: 100/500/2000/5000 cells, 1 s each.
    srun --mpi=pmix "$PY" neuron/tc_mpi.py --bench --gh-tc "$GH" --gsk-re "$GSK"
else
    srun --mpi=pmix "$PY" neuron/tc_mpi.py \
        --n-tc "$NTC" --n-re "$NRE" --tstop "$TSTOP" \
        --gh-tc "$GH" --gsk-re "$GSK" \
        --out "out/${TAG}.npz"
fi

echo "== done; results in out/ =="
