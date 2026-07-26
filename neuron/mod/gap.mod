: Electrical (gap-junction) coupling between reticular (TRN) cells via
: connexin-36 (Fernandez & Luthi 2020, sect. V.C.1: gap junctions synchronise
: spindle-like rhythms within a TRN sector). A POINTER reads the partner cell's
: membrane potential; the coupling current is ohmic. Bidirectional coupling is
: made by placing one Gap on each cell, each pointing at the other's v. Suited
: to the weak coupling of TRN (strong coupling would want a simultaneous solve).
NEURON {
    POINT_PROCESS Gap
    NONSPECIFIC_CURRENT i
    RANGE g, i
    POINTER vgap
}
UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}
PARAMETER { g = 0 (uS) }
ASSIGNED { v (mV) vgap (mV) i (nA) }
BREAKPOINT { i = g * (v - vgap) }
