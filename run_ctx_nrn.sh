#!/bin/bash
#SBATCH --job-name=audcc_ctx_nrn
#SBATCH --output=slurmout_ctx_%j.txt
#SBATCH --error=slurmerr_ctx_%j.txt
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --qos=gp_debug
#
# MareNostrum 5 -- MPI-parallel NEURON CORTICOTHALAMIC network (thalamic
# TC<->RE loop + a full L4/L2-3/L5/L6 cortical column: hh2 + it/it2 + SK2 +
# Ca-dependent I_h in thalamus, hh2 + ical + inap + SK2/Ca2+ adaptation in
# cortex). Companion to run_nrn.sh (thalamus only) -- use this one once the
# thalamus-only benchmark/production runs are validated and you want the
# closed cortex->thalamus->cortex loop at the reference model's bio-plausible
# size (config/network_auditory_mn5_5k.yaml: 3050 declared cells).
#
#   sbatch --export=ALL,BENCH=1 run_ctx_nrn.sh
#       -> scaling benchmark: 0.1x/0.5x/1x/1.65x the default 3050-cell sizes
#          (1.65x ~= 5010, matching the NEST reference model's actual size
#          after its internal subtype splits)
#
#   sbatch --export=ALL,SCALE=1.0,TSTOP=200000 run_ctx_nrn.sh
#       -> production at the reference size (SCALE=1.0 -> 3050 cells;
#          SCALE=1.65 -> ~5010)
#
# RUN THE BENCHMARK FIRST -- same reasoning as run_nrn.sh: this tells you
# whether the requested scale is affordable in the time/rank budget before a
# 200 s production run finds out the hard way. Also run run_nrn.sh's own
# --bench once if you have not already: the thalamus-only benchmark isolates
# whether any slowdown comes from the cortex or from the thalamic loop.

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$WORKDIR"
[ -f neuron/ctx_thalamus_mpi.py ] || { echo "ERROR: run sbatch from the repo root (no neuron/ctx_thalamus_mpi.py in $WORKDIR)"; exit 1; }

# --- environment -------------------------------------------------------
# Same convention as run_nrn.sh: this script loads NO modules. Load
# `module load miniforge` (or whatever provides python + NEURON on your
# MN5 project) BEFORE sbatch, and leave this script alone -- see run_nrn.sh's
# comments for the module-order pitfalls that were hit and fixed there.

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

# --- compile the NMODL mechanisms (whole mod/ dir: it, it2, hh2, cad, sk2,
#     ihca, ical, inap, gap, gapmpi) ------------------------------------
cd neuron
if [ ! -f x86_64/libnrnmech.so ]; then
    echo "== compiling NMODL mechanisms =="
    command -v nrnivmodl >/dev/null || { echo "ERROR: nrnivmodl not on PATH -- is the NEURON module loaded?"; exit 1; }
    nrnivmodl mod
fi
cd "$WORKDIR"

SCALE="${SCALE:-1.0}"
TSTOP="${TSTOP:-200000}"
CONV="${CONV:-100}"
TAG="${TAG:-ctx_nrn_${SLURM_JOB_ID:-local}}"
mkdir -p out

echo "== NEURON MPI (corticothalamic): ${SLURM_NTASKS:-1} ranks, scale=${SCALE} =="

if [ -n "${SWEEP_G_RE_RE:-}" ]; then
    # SWEEP_G_RE_RE: comma-separated g_re_re values, e.g. "0.005,0.01,0.02,0.05,0.1"
    # -- runs each at full production --scale (SCALE, default 1.0) in one job,
    # reporting RE burst SHAPE (frac_burst, mean_burst_size, event_hz), not
    # just spike counts. SWEEP_TSTOP (optional, default 20000 ms) keeps each
    # point cheap; this is NOT the 200s production run.
    SWEEP_ARGS=()
    [ -n "${SWEEP_TSTOP:-}" ] && SWEEP_ARGS+=(--sweep-tstop "$SWEEP_TSTOP")
    srun --mpi=pmix "$PY" neuron/ctx_thalamus_mpi.py \
        --sweep-g-re-re "$SWEEP_G_RE_RE" --scale "$SCALE" --conv "$CONV" \
        "${SWEEP_ARGS[@]}"
elif [ "${BENCH:-0}" = "1" ]; then
    # BENCH_SCALES (optional): comma-separated --scale values, e.g. "0.02,0.05"
    # for a fast small-rank correctness check. Unset -> ctx_thalamus_mpi.py's
    # own default (0.1,0.5,1.0,1.65).
    BENCH_ARGS=()
    [ -n "${BENCH_SCALES:-}" ] && BENCH_ARGS+=(--bench-scales "$BENCH_SCALES")
    [ -n "${BENCH_MS:-}" ] && BENCH_ARGS+=(--bench-ms "$BENCH_MS")
    srun --mpi=pmix "$PY" neuron/ctx_thalamus_mpi.py --bench --conv "$CONV" "${BENCH_ARGS[@]}"
else
    # G_RE_RE/G_RE_RE_SD (optional): mean/SD of the distributed RE<->RE
    # lateral-inhibition weight. Unset -> ctx_thalamus_mpi.py's own defaults
    # (0.006/0.002, from the sweep in job 45081779).
    PROD_ARGS=()
    [ -n "${G_RE_RE:-}" ] && PROD_ARGS+=(--g-re-re "$G_RE_RE")
    [ -n "${G_RE_RE_SD:-}" ] && PROD_ARGS+=(--g-re-re-sd "$G_RE_RE_SD")
    srun --mpi=pmix "$PY" neuron/ctx_thalamus_mpi.py \
        --scale "$SCALE" --tstop "$TSTOP" --conv "$CONV" \
        --out "out/${TAG}.npz" "${PROD_ARGS[@]}"
fi

echo "== done; results in out/ =="
