"""
Brain-state presets: NREM sleep vs wake as neuromodulation of a K+ leak
(PLAN-auditory-input.md, design D4).

Every cell's passive leak is split into  pas + kleak  (mod/kleak.mod):

    g_tot * (v - e_old)  ==  g_pas * (v - e_open) + g_kl * (v - E_KL)

with the SAME total conductance g_tot and the SAME combined reversal e_old
(the cell's current, sleep-state e_pas). So the "nrem" state is the existing
model, unchanged up to floating-point rounding. e_open is the reversal the
cell relaxes to when neuromodulators close the K+ leak completely: the wake
direction is `kl_scale` < 1, which depolarises the cell and raises its input
resistance (McCormick 1992; Bazhenov et al. 2002).

Stage A only defines and applies the leak split. The other wake knobs in D4
(cortical SK2 / I_NaP scaling, L5 gap scaling, background noise) and the
calibration of the "wake" values below belong to Stage B, so the "wake"
preset here is a placeholder: it closes the K+ leak and nothing else.
"""

E_KL = -95.0   # K+ leak reversal (mV)

# Reversal with the K+ leak fully closed ("open" = neuromodulated), per cell
# group. Starting points for Stage B, not calibrated values: awake TC cells
# rest around -60 to -65 mV (tonic relay mode, I_T inactivated); RE somewhat
# lower; cortical cells a few mV above their current -70 / -68 mV.
E_OPEN = {"tc": -62.0, "re": -68.0, "cx_e": -64.0, "cx_i": -64.0}

# kl_scale multiplies each group's K+ leak conductance after the split.
STATES = {
    "nrem": {"kl_scale": {"tc": 1.0, "re": 1.0, "cx_e": 1.0, "cx_i": 1.0}},
    # PLACEHOLDER until Stage B calibrates it (see module docstring). As is,
    # it makes an isolated TC cell fire at ~35 Hz with no input
    # (aud_checks.py --only tc): the TCCell fires tonically from about -67 mV,
    # so E_OPEN["tc"] = -62 with the K+ leak fully closed is too depolarising.
    "wake": {"kl_scale": {"tc": 0.0, "re": 0.0, "cx_e": 0.0, "cx_i": 0.0}},
}


def group_of(pop):
    """Map a ctx_thalamus_mpi population name to a leak group."""
    if pop in ("tc", "re"):
        return pop
    return "cx_i" if pop.endswith("i") else "cx_e"


def split_leak(sec, e_open, kl_scale=1.0, ekl=E_KL):
    """Replace sec's pas leak by pas + kleak with the same total conductance
    and combined reversal, then scale the K+ part by kl_scale.

    Returns the K+ fraction f of the total leak (0 if e_pas is already at or
    above e_open, in which case no K+ leak is added)."""
    e_old = sec(0.5).e_pas
    g_tot = sec(0.5).g_pas
    f = (e_old - e_open) / (ekl - e_open)
    f = min(max(f, 0.0), 0.99)
    sec.insert("kleak")
    sec.ekl_kleak = ekl
    sec.g_kleak = f * g_tot * kl_scale
    sec.g_pas = (1.0 - f) * g_tot
    sec.e_pas = e_open
    return f


def apply_state(cell, pop, state):
    """Split the leak in every section of `cell` that has pas."""
    grp = group_of(pop)
    scale = STATES[state]["kl_scale"][grp]
    for name in ("soma", "dend"):
        sec = getattr(cell, name, None)
        if sec is not None and sec.has_membrane("pas"):
            split_leak(sec, E_OPEN[grp], scale)
