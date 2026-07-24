: Submembrane Ca2+ pool (Destexhe et al. 1993): accumulates from the T-current,
: decays to baseline, provides cai to SK2. Concentrations are in NEURON's native
: mM; molar-family unit names are avoided (NEURON 9 does not define them).
NEURON {
    SUFFIX cad
    USEION ca READ ica, cai WRITE cai
    RANGE depth, taur, cainf
}
UNITS {
    (mA) = (milliamp)
    FARADAY = (faraday) (coulomb)
}
PARAMETER {
    depth = 0.1
    taur  = 80  (ms)
    cainf = 5e-5
}
ASSIGNED { ica (mA/cm2) drive }
STATE { cai }
BREAKPOINT { SOLVE state METHOD cnexp }
DERIVATIVE state {
    drive = -(10000) * ica / (2 * FARADAY * depth)
    if (drive <= 0) { drive = 0 }
    cai' = drive + (cainf - cai)/taur
}
INITIAL { cai = cainf }
