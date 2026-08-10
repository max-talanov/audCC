# Uploading to MareNostrum 5

Everything needed for a NEST run is **~196 KB** — the model is pure Python, so
only source and configs go up. Outputs are generated on MN5.

## Files to upload

| file | why |
|------|-----|
| `run.sh` | the SLURM job script (`sbatch run.sh`) |
| `tc_network.py` | the model — builds/wires/drives the column |
| `tc_run.py` | driver; **self-validates** (non-zero exit if a rhythm is missing) |
| `requirements.txt` | dependency list (NEST is installed separately) |
| `config/*.yaml` | all six configs, incl. `network_auditory_mn5_5k.yaml` (~5000 neurons, HH + SK2) |
| `tc_validate.py` | scores vs Fernandez & Lüthi (2020) criteria — needed if `VALIDATE=1` |
| `tc_spindle_figures.py` | **required by `tc_validate.py`** (imports `tc_episodes`) |
| `tc_analyze.py` | per-layer spindle metrics + propagation lag — needed if `ANALYZE=1` |
| `tc_present.py` | presentation figure — needed if `PRESENT=1` |
| `tc_io.py` | **required** — writes the `.h5` raw result (`SAVE_H5=1`, the default) |
| `tc_architecture.py` | architecture schematic (optional; no NEST needed) |
| `__init__.py` | package marker |

## Do NOT upload

| path | size | why |
|------|------|-----|
| `.venv-neuron/` | 495 MB | local NEURON venv; rebuild on MN5 if ever needed |
| `out/` | 7.1 MB | generated outputs |
| `neuron/` | 748 KB | the NEURON port — **not used** for NEST runs on MN5 |
| `docs/`, `*.pdf` | — | documentation / reference paper, not needed at runtime |
| `.git/` | — | use `git clone` instead if you want history |

## scp

One-liner from the repo root (creates the remote dir first):

```bash
ssh USER@glogin1.bsc.es 'mkdir -p ~/audCC/config'
scp run.sh requirements.txt tc_network.py tc_run.py tc_io.py tc_validate.py \
    tc_spindle_figures.py tc_analyze.py tc_present.py tc_architecture.py \
    __init__.py  USER@glogin1.bsc.es:~/audCC/
scp config/*.yaml  USER@glogin1.bsc.es:~/audCC/config/
```

Or as a single tarball (fewer round trips):

```bash
tar czf audcc.tgz run.sh requirements.txt tc_*.py __init__.py config/*.yaml
scp audcc.tgz USER@glogin1.bsc.es:~/
ssh USER@glogin1.bsc.es 'mkdir -p ~/audCC && tar xzf audcc.tgz -C ~/audCC'
```

Or, since the repo is on GitHub, skip scp entirely:

```bash
ssh USER@glogin1.bsc.es 'git clone https://github.com/max-talanov/audCC.git'
```

## On MN5, before the first run

1. **Provide NEST 3.x** — it is *not* pip-installable. Load a module or activate
   an environment that supplies it, and uncomment the matching lines in the
   `# ---- environment ----` block of `run.sh`.
2. Install the plain Python deps into that environment:
   `pip install numpy scipy matplotlib pyyaml`
   (`requirements.txt` also lists `neuron`/`tvb-*`, which the NEST runs do not need.)
3. Set `--account` / `--qos` / `--partition` in `run.sh` to your project.
4. Smoke test before the production run:
   ```bash
   sbatch --export=ALL,CONFIG=config/network_auditory_mn5_5k.yaml,TSTOP=3000 run.sh
   ```
   ~1 min, 5010 neurons. Then the real run with `TSTOP=200000,VALIDATE=1,ANALYZE=1`.

## Getting results back

`SAVE_H5=1` (the default) writes the raw spikes/traces to `out/<TAG>.h5`. Pull it
back and regenerate **every** figure locally — no NEST, no re-simulation:

```bash
scp USER@glogin1.bsc.es:'~/audCC/out/*.h5' .
python3 tc_plot_h5.py mn5_5k_200s.h5 --outdir out     # figures + validation + analysis
python3 tc_plot_h5.py mn5_5k_200s.h5 --info           # just summarise the file
```

`tc_plot_h5.py` runs locally only and needs just numpy / scipy / matplotlib / h5py.
