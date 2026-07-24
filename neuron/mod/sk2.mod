: SK2 (small-conductance Ca2+-activated K+), gated by the PRIVATE SK Ca pool
: ("sk" ion from cad) via a Hill relation. Voltage-independent. Terminates the
: low-threshold Ca2+ spike after a few fast spikes -- the burst after-
: hyperpolarisation that keeps thalamic bursts SHORT (Fernandez & Luthi V.A.1).
NEURON {
    SUFFIX sk2
    USEION sk READ ski VALENCE 2
    USEION k READ ek WRITE ik
    RANGE gkbar, g, kd, hill, tauz
}
UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S)  = (siemens)
}
PARAMETER {
    gkbar = 0.002 (S/cm2)
    kd    = 0.5
    hill  = 4
    tauz  = 12 (ms)
}
ASSIGNED { v (mV) ek (mV) ski ik (mA/cm2) g (S/cm2) zinf }
STATE { z }
BREAKPOINT {
    SOLVE state METHOD cnexp
    g = gkbar * z
    ik = g * (v - ek)
}
DERIVATIVE state {
    zinf = 1/(1 + (kd/ski)^hill)
    z' = (zinf - z)/tauz
}
INITIAL {
    zinf = 1/(1 + (kd/ski)^hill)
    z = zinf
}
