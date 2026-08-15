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

module purge
module load python || true
# module load openmpi   # uncomment if NEURON's MPI needs the system stack

PY="${PY:-./.venv-neuron/bin/python}"
[ -x "$PY" ] || { echo "ERROR: no python at $PY -- create .venv-neuron and pip install neuron"; exit 1; }

# Compile the NMODL mechanisms next to the sources (itd, hh2, cad, sk2, ihca, gap)
cd neuron
if [ ! -f x86_64/libnrnmech.so ]; then
    echo "== compiling NMODL mechanisms =="
    "../$(dirname "$PY")/nrnivmodl" mod || ../.venv-neuron/bin/nrnivmodl mod
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
