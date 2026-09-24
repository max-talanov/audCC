TITLE K+ leak current -- the neuromodulatory sleep/wake knob

COMMENT
Voltage-independent K+ leak, i = g * (v - ekl).

This is the conductance ACh / NE / 5-HT / histamine close in wakefulness
(McCormick 1992), depolarising thalamic cells out of burst mode and raising
input resistance. Bazhenov et al. (2002) and Hill & Tononi (2005, g_KL) use it
as the model's sleep-depth knob; PLAN-auditory-input.md D4 does the same.

A NONSPECIFIC_CURRENT with its own reversal (ekl), not USEION k: hh2 and sk2
READ ek, and this current must not change their driving force.

ctx_thalamus_mpi.py splits each section's existing pas leak into pas + kleak
with the same total conductance and the same combined reversal, so the NREM
state is unchanged; scaling g down is the wake direction.
ENDCOMMENT

NEURON {
    SUFFIX kleak
    NONSPECIFIC_CURRENT i
    RANGE g, ekl
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
}

PARAMETER {
    g   = 0     (S/cm2)
    ekl = -95   (mV)
}

ASSIGNED {
    v (mV)
    i (mA/cm2)
}

BREAKPOINT {
    i = g * (v - ekl)
}
