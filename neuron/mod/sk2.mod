: SK2 (small-conductance Ca2+-activated K+). Voltage-INDEPENDENT, gated by
: intracellular Ca2+ via a Hill relation. Fernandez & Luthi 2020 (sect. V.A.1):
: SK2 produces the burst after-hyperpolarisation that keeps thalamic bursts
: SHORT and repeatable -- i.e. it terminates the low-threshold Ca2+ spike after
: 2-6 fast spikes instead of the ~250 ms runaway plateau seen without it.
NEURON {
    SUFFIX sk2
    USEION k READ ek WRITE ik
    USEION ca READ cai
    RANGE gkbar, g, kd
}
UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S)  = (siemens)
}
PARAMETER {
    gkbar = 0.002 (S/cm2)
    kd    = 0.0007        : half-activation Ca2+ (mM, native)
    hill  = 4
    tauz  = 12 (ms)
}
ASSIGNED { v (mV) ek (mV) cai ik (mA/cm2) g (S/cm2) zinf }
STATE { z }
BREAKPOINT {
    SOLVE state METHOD cnexp
    g = gkbar * z
    ik = g * (v - ek)
}
DERIVATIVE state {
    zinf = 1/(1 + (kd/cai)^hill)
    z' = (zinf - z)/tauz
}
INITIAL {
    zinf = 1/(1 + (kd/cai)^hill)
    z = zinf
}
