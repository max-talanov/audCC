#!/bin/bash -l
#SBATCH --job-name=AUDCC_NEST
#SBATCH --output=audcc_%A_%a.slurmout
#SBATCH --error=audcc_%A_%a.slurmerr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=02:00:00
#SBATCH --partition=gp_bsccs

# Auditory thalamo-cortical sleep model (slow waves + sleep spindles) on MN5.
# One task, many NEST threads (shared-memory); NEST scales ~linearly in neurons.
#
# ---- typical runs -----------------------------------------------------------
#
# 5k neurons, Hodgkin-Huxley WITH SK2, 3 s smoke test (~1 min)
#  sbatch --export=ALL,CONFIG=config/network_auditory_mn5_5k.yaml,TSTOP=3000 run.sh
#
# 5k neurons, HH + SK2, 200 s production run (spindle statistics; ~35 min)
#  sbatch --export=ALL,CONFIG=config/network_auditory_mn5_5k.yaml,TSTOP=200000,VALIDATE=1,ANALYZE=1 run.sh
#
# The 10/10-validated AdEx spindle model (needs TRIGGER=1)
#  sbatch --export=ALL,CONFIG=config/network_auditory_adex.yaml,TSTOP=200000,TRIGGER=1,VALIDATE=1,PRESENT=1 run.sh
#
# Classic iaf column (fast, full cortical activity, imposed-drive spindles)
#  sbatch --export=ALL,CONFIG=config/network_auditory_mn5.yaml,TSTOP=30000,ANALYZE=1 run.sh
#
# Force a neuron model regardless of what the config says
#  sbatch --export=ALL,CONFIG=config/network_auditory_mn5.yaml,MODEL=ht_neuron run.sh
#
# ---- notes ------------------------------------------------------------------
#  * Adjust --account/--qos/--partition to your MN5 project.
#  * MN5 gpp nodes have 112 cores; --cpus-per-task can go up to that. Measured
#    thread scaling on this model tapered past ~8 threads on a laptop, so more
#    cores help but sub-linearly -- do not assume 112x.
#  * tc_run.py SELF-VALIDATES: it exits non-zero if the ~1 Hz slow wave or the
#    spindle is missing, so a failed job is a real scientific failure, not just
#    a crash.
#  * TSTOP >= 150000 (150 s) is needed for spindle density / inter-spindle
#    interval / the ~0.02 Hz infraslow clustering. 3 s is a smoke test only.

CONFIG=${CONFIG:-config/network_auditory_mn5_5k.yaml}
TSTOP=${TSTOP:-}                # ms; empty = use the value in the YAML
MODEL=${MODEL:-}                # empty = use simulation.neuron_model from YAML
SEED=${SEED:-1}
TRIGGER=${TRIGGER:-0}           # 1 = external modulatory spindle trigger (HH/AdEx)
VALIDATE=${VALIDATE:-0}         # 1 = score vs Fernandez & Luthi 2020 criteria
ANALYZE=${ANALYZE:-0}           # 1 = per-layer spindle metrics + propagation lag
PRESENT=${PRESENT:-0}           # 1 = presentation figure (slow waves + spindles)
NO_PLOT=${NO_PLOT:-0}           # 1 = skip figures (faster)
OUTDIR=${OUTDIR:-out}
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

export LANG=${LANG:-C.UTF-8}
export LC_ALL=${LC_ALL:-C.UTF-8}
export PYTHONIOENCODING=utf-8
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg           # headless plotting on compute nodes

unset OMP_NUM_THREADS
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# ---- environment ------------------------------------------------------------
# Load a Python/NEST module or activate a venv that provides NEST 3.x here.
# NEST is NOT pip-installable (see requirements.txt).
# module load python/3.12
# module load nest/3.9
# source "$HOME/audcc-venv/bin/activate"

echo "[Slurm] job=$SLURM_JOB_ID  ntasks=$SLURM_NTASKS  cpus-per-task=$SLURM_CPUS_PER_TASK"
echo "[Slurm] config=${CONFIG}  tstop=${TSTOP:-<from yaml>}  model=${MODEL:-<from yaml>}  seed=${SEED}"
echo "[Slurm] trigger=${TRIGGER}  validate=${VALIDATE}  analyze=${ANALYZE}  present=${PRESENT}  no_plot=${NO_PLOT}"

python3 - <<'PY'
import nest
ks = nest.GetKernelStatus()
mpi = ks.get("mpi_num_processes", ks.get("num_processes", ks.get("total_num_processes", 1)))
thr = ks.get("local_num_threads", ks.get("num_threads", ks.get("threads", 1)))
print(f"nest {nest.__version__}  mpi_procs={mpi}  local_threads={thr}")
PY

mkdir -p "$OUTDIR" logs

# ---- output tag from the active options --------------------------------------
CFG_TAG=$(basename "$CONFIG" .yaml | sed 's/^network_auditory_//')
TAG="${CFG_TAG}"
[ -n "$MODEL" ]     && TAG="${TAG}_${MODEL}"
[ "$TRIGGER" = "1" ] && TAG="${TAG}_trig"
[ -n "$TSTOP" ]     && TAG="${TAG}_$((TSTOP/1000))s"
echo "[Slurm] tag → ${TAG}   (figures/reports land in ${OUTDIR}/)"

# ---- optional flags ----------------------------------------------------------
FLAGS=""
[ -n "$TSTOP" ]      && FLAGS="$FLAGS --tstop $TSTOP"
[ -n "$MODEL" ]      && FLAGS="$FLAGS --neuron-model $MODEL"
[ "$TRIGGER" = "1" ] && FLAGS="$FLAGS --spindle-trigger"
[ "$NO_PLOT" = "1" ] && FLAGS="$FLAGS --no-plot"

# ---- main run (self-validating: non-zero exit if a rhythm is missing) --------
srun --cpu-bind=cores \
  python3 -u tc_run.py \
    --config  "$CONFIG" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --outdir  "$OUTDIR" \
    --tag     "$TAG" \
    --seed    "$SEED" \
    $FLAGS
RC=$?
echo "[Slurm] tc_run.py exit=$RC"

# ---- optional post-processing ------------------------------------------------
if [ "$VALIDATE" = "1" ]; then
  echo "[Slurm] validating against Fernandez & Luthi 2020 criteria..."
  VFLAGS=""
  [ -n "$TSTOP" ] && VFLAGS="$VFLAGS --tstop $TSTOP"
  [ -n "$MODEL" ] && VFLAGS="$VFLAGS --neuron-model $MODEL"
  [ "$TRIGGER" != "1" ] && VFLAGS="$VFLAGS --no-trigger"
  srun --cpu-bind=cores \
    python3 -u tc_validate.py --config "$CONFIG" --seed "$SEED" \
      $VFLAGS --no-assert
fi

if [ "$ANALYZE" = "1" ]; then
  echo "[Slurm] per-layer spindle analysis..."
  srun --cpu-bind=cores \
    python3 -u tc_analyze.py --config "$CONFIG" --seed "$SEED" \
      ${TSTOP:+--tstop "$TSTOP"} --outdir "$OUTDIR"
fi

if [ "$PRESENT" = "1" ]; then
  echo "[Slurm] presentation figure..."
  srun --cpu-bind=cores \
    python3 -u tc_present.py --config "$CONFIG" --seed "$SEED" \
      ${TSTOP:+--tstop "$TSTOP"} --outdir "$OUTDIR"
fi

echo "[Slurm] done. outputs in ${OUTDIR}/"
exit $RC
