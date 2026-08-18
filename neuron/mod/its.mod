: Low-threshold T-type Ca2+ current for RETICULAR (TRN/nRT) cells -- Ca_v3.3.
: Destexhe, Contreras, Steriade, Sejnowski & Huguenard (1996), ITs.mod;
: kinetics from Huguenard & Prince (1992) reticular recordings.
:
: The relay cell's current (itd.mod) is Ca_v3.1. The reticular nucleus expresses
: Ca_v3.3, which activates at a slightly more depolarised potential and -- the
: point here -- has DIFFERENT inactivation kinetics. Using one mechanism for
: both cell types (as this model did) is biologically wrong and forces both
: populations onto the same cycle period. TRN is the spindle pacemaker
: (Fernandez & Luthi 2020 sect. V), so its kinetics set the loop rhythm.
:
: GHK driving force, temperature-scaled taus, same conventions as itd.mod.

NEURON {
    SUFFIX its
    USEION ca READ cai, cao WRITE ica
    RANGE pcabar, m_inf, tau_m, h_inf, tau_h, shift, ica
}

UNITS {
    (molar) = (1/liter)
    (mV)    = (millivolt)
    (mA)    = (milliamp)
    (mM)    = (millimolar)
    FARADAY = (faraday) (coulomb)
    R       = (k-mole)  (joule/degC)
}

PARAMETER {
    pcabar = 2.5e-4 (cm/s)
    shift  = 2      (mV)
    celsius (degC)
}

ASSIGNED {
    v (mV)  cai (mM)  cao (mM)  ica (mA/cm2)
    m_inf  tau_m (ms)  h_inf  tau_h (ms)  phi_m  phi_h
}

STATE { m h }

BREAKPOINT {
    SOLVE castate METHOD cnexp
    ica = pcabar * m*m*h * ghk(v, cai, cao)
}

DERIVATIVE castate {
    evaluate_fct(v)
    m' = (m_inf - m) / tau_m
    h' = (h_inf - h) / tau_h
}

INITIAL {
    phi_m = 5.0 ^ ((celsius - 24)/10)
    phi_h = 3.0 ^ ((celsius - 24)/10)
    evaluate_fct(v)
    m = m_inf
    h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL vs
    vs = v + shift
    m_inf = 1.0 / (1 + exp(-(vs + 52)/7.4))
    h_inf = 1.0 / (1 + exp((vs + 80)/5))
    tau_m = (3 + 1.0/(exp((vs + 27)/10) + exp(-(vs + 102)/15))) / phi_m
    tau_h = (85 + 1.0/(exp((vs + 48)/4) + exp(-(vs + 407)/50))) / phi_h
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) { LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    eco = co*efun(z)
    eci = ci*efun(-z)
    ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
    if (fabs(z) < 1e-4) { efun = 1 - z/2 }
    else { efun = z/(exp(z) - 1) }
}
