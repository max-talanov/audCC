"""Noise floor for the regression gate: legacy network with a 1e-9 mV
perturbation of e_pas (PERT_MODE=all: every cell; PERT_MODE=one: gid 0 only),
i.e. the size of the rounding difference the pas + kleak split introduces.

Run from neuron/ (mechanisms compiled there), same flags as ctx_thalamus_mpi.py:
    PERT_MODE=all mpiexec -n 4 ../.venv-neuron/bin/python \
        ../res/2026-09-24/stageA_perturbed.py --scale 0.1 --tstop 20000 \
        --g-i-e-l5 0.02 --g-l5-rec 0.013 --tau2-l5-rec 200 --taur-l5e-rs 120 \
        --l5-rec-mech nmda --out legacy_pert_all.npz
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "neuron"))
import ctx_thalamus_mpi as M
mode = os.environ["PERT_MODE"]
orig = M.ParallelCorticoThalamicNet._make_cell
def make(self, pop, gid):
    c = orig(self, pop, gid)
    if mode == "all" or gid == 0:
        c.soma.e_pas += 1e-9
    return c
M.ParallelCorticoThalamicNet._make_cell = make
M.main()
