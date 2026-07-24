: Low-threshold T-type Ca2+ current (Ca_v3.1 / Ca_v3.3), Destexhe-style.
: The current behind thalamic rebound bursting (Fernandez & Luthi 2020, sect. V):
: de-inactivates when the cell is hyperpolarised, so release from inhibition
: produces a low-threshold Ca2+ spike that fires a burst of Na+ spikes.
NEURON {
    SUFFIX cav3
    USEION ca READ cai, cao WRITE ica
    RANGE gmax, g, ica
}
UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S)  = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}
PARAMETER {
    gmax = 0.002 (S/cm2)
    celsius (degC)
}
ASSIGNED {
    v (mV)
    cai (mM)
    cao (mM)
    ica (mA/cm2)
    g (S/cm2)
    minf
    hinf
    mtau (ms)
    htau (ms)
}
STATE { m h }
BREAKPOINT {
    SOLVE states METHOD cnexp
    : ghk-free ohmic form is adequate here; drive by (v - eca-ish) via Nernst
    g = gmax * m*m*h
    ica = g * (v - eca())
}
FUNCTION eca() (mV) {
    eca = (1000)*(R*(celsius+273.15))/(2*FARADAY)*log(cao/cai)
}
INITIAL {
    rates(v)
    m = minf
    h = hinf
}
DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
    h' = (hinf - h)/htau
}
PROCEDURE rates(v(mV)) {
    : activation (fast) and inactivation after Huguenard & McCormick 1992.
    : htau is TWO-BRANCH: hyperpolarised de-inactivation is a few hundred ms
    : (recovery), NOT thousands -- a single-branch form explodes below -80 mV
    : and prevents rebound within a physiological hyperpolarisation.
    minf = 1/(1 + exp(-(v + 57)/6.2))
    hinf = 1/(1 + exp((v + 81)/4.0))
    mtau = 0.612 + 1/(exp(-(v + 132)/16.7) + exp((v + 16.8)/18.2))
    if (v < -80) {
        htau = exp((v + 467)/66.6)
    } else {
        htau = 28 + exp(-(v + 22)/10.5)
    }
}
