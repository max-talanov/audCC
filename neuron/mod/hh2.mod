: Fast Na+ and delayed-rectifier K+ (Traub & Miles 1991), as used in Destexhe's
: thalamocortical models. Unlike NEURON's built-in `hh`, these kinetics recover
: quickly enough to fire a REPETITIVE burst on a sustained depolarisation (the
: T-current Ca2+ plateau) instead of going into depolarisation block after one
: spike. `vtraub` shifts threshold to the resting range.
NEURON {
    SUFFIX hh2
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    RANGE gnabar, gkbar, vtraub
}
UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}
PARAMETER {
    gnabar = 0.05 (mho/cm2)
    gkbar  = 0.005 (mho/cm2)
    vtraub = -63 (mV)
    celsius (degC)
}
ASSIGNED {
    v (mV)
    ena (mV)
    ek (mV)
    ina (mA/cm2)
    ik (mA/cm2)
    m_inf h_inf n_inf
    tau_m tau_h tau_n (ms)
    tadj
}
STATE { m h n }
BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gnabar * m*m*m*h * (v - ena)
    ik  = gkbar  * n*n*n*n * (v - ek)
}
DERIVATIVE states {
    evaluate_fct(v)
    m' = (m_inf - m)/tau_m
    h' = (h_inf - h)/tau_h
    n' = (n_inf - n)/tau_n
}
INITIAL {
    tadj = 3.0 ^ ((celsius - 36)/10)
    evaluate_fct(v)
    m = m_inf
    h = h_inf
    n = n_inf
}
PROCEDURE evaluate_fct(v(mV)) { LOCAL a, b, v2
    v2 = v - vtraub
    a = 0.32 * vtrap(13 - v2, 4)
    b = 0.28 * vtrap(v2 - 40, 5)
    tau_m = 1 / (a + b) / tadj
    m_inf = a / (a + b)

    a = 0.128 * exp((17 - v2)/18)
    b = 4 / (1 + exp((40 - v2)/5))
    tau_h = 1 / (a + b) / tadj
    h_inf = a / (a + b)

    a = 0.032 * vtrap(15 - v2, 5)
    b = 0.5 * exp((10 - v2)/40)
    tau_n = 1 / (a + b) / tadj
    n_inf = a / (a + b)
}
FUNCTION vtrap(x, y) {
    if (fabs(x/y) < 1e-6) { vtrap = y * (1 - x/y/2) }
    else { vtrap = x / (exp(x/y) - 1) }
}
