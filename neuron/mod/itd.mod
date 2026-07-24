: Low-threshold T-type Ca2+ current, Destexhe, Bal, McCormick & Sejnowski 1996
: (J Neurophysiol 76:2049) -- the classic thalamic-relay IT used across the
: Destexhe ModelDB thalamocortical models. GHK driving force + temperature-
: scaled (phi) kinetics, so the low-threshold Ca2+ spike is a BRIEF (~60-100 ms)
: event carrying 2-6 fast spikes -- unlike a fixed-ohmic form with unscaled taus,
: which gives an over-long plateau. Molar unit names are omitted (NEURON 9 does
: not define them); cai/cao use NEURON's native mM.
NEURON {
    SUFFIX itd
    USEION ca READ cai, cao WRITE ica
    RANGE pcabar, shift, m, h
    GLOBAL qm, qh
}
UNITS {
    (mV)  = (millivolt)
    (mA)  = (milliamp)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}
PARAMETER {
    pcabar = 2.5e-4      : permeability (cm/s scale)
    shift  = 2 (mV)
    qm = 3.55
    qh = 3.0
    celsius (degC)
}
ASSIGNED {
    v (mV) cai cao ica (mA/cm2)
    phim phih minf hinf taum (ms) tauh (ms)
}
STATE { m h }
BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = pcabar * m*m*h * ghk(v, cai, cao)
}
DERIVATIVE states {
    rates(v)
    m' = (minf - m)/taum
    h' = (hinf - h)/tauh
}
INITIAL {
    phim = qm^((celsius-24)/10)
    phih = qh^((celsius-24)/10)
    rates(v)
    m = minf
    h = hinf
}
PROCEDURE rates(v(mV)) { LOCAL vs
    vs = v + shift
    minf = 1/(1 + exp(-(vs+57)/6.2))
    hinf = 1/(1 + exp((vs+81)/4))
    taum = (0.612 + 1/(exp(-(vs+132)/16.7) + exp((vs+16.8)/18.2))) / phim
    if (vs < -80) {
        tauh = exp((vs+467)/66.6) / phih
    } else {
        tauh = (28 + exp(-(vs+22)/10.5)) / phih
    }
}
FUNCTION ghk(v(mV), ci, co) {
    LOCAL z
    z = (2e-3)*FARADAY*v/(R*(celsius+273.15))
    ghk = (2e-3)*FARADAY*(ci*efun(-z) - co*efun(z))
}
FUNCTION efun(z) {
    if (fabs(z) < 1e-4) { efun = 1 - z/2 }
    else { efun = z/(exp(z) - 1) }
}
