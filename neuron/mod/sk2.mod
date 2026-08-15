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
    : Half-activating Ca2+ concentration, in NEURON's native mM. SK channels
    : have K_d ~= 0.3-0.7 uM (Hirschberg et al. 1998; Xia et al. 1998), i.e.
    : 3e-4 - 7e-4 mM. NOTE: this was 0.5 -- a uM value left in a mM field, so
    : 1000x too high. With the submembrane pool peaking near 3.7e-2 mM and
    : hill=4, that held activation at ~3e-5 of gkbar: the channel never opened,
    : which is why earlier sweeps of gkbar (0 -> 0.03) changed nothing.
    kd    = 0.0005
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
