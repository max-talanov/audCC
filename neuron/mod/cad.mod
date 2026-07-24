: Submembrane Ca2+ pool for SK2, in a PRIVATE ion species ("sk") so it never
: feeds back into the T-current's Ca2+ reversal. It READs the real Ca current
: (ica, from cav3) as the Ca source but WRITEs a separate concentration that
: only sk2 reads. This decouples the SK microdomain from the bulk Ca that sets
: the T-current driving force (Fernandez & Luthi sect. V.A.1; the microdomain
: and bulk Ca pools are physiologically distinct).
NEURON {
    SUFFIX cad
    USEION ca READ ica
    USEION sk WRITE ski VALENCE 2
    RANGE depth, taur, skinf
}
UNITS {
    (mA) = (milliamp)
    FARADAY = (faraday) (coulomb)
}
PARAMETER {
    depth = 0.1
    taur  = 80 (ms)
    skinf = 5e-5
}
ASSIGNED { ica (mA/cm2) drive }
STATE { ski }
BREAKPOINT { SOLVE state METHOD cnexp }
DERIVATIVE state {
    drive = -(10000) * ica / (2 * FARADAY * depth)
    if (drive <= 0) { drive = 0 }
    ski' = drive + (skinf - ski)/taur
}
INITIAL { ski = skinf }
