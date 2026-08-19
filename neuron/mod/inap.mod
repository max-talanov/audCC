TITLE Persistent Na+ current (cortical layer-5 intrinsic bursting)
:
: Instantaneous m_inf-gated persistent Na+ current (Compte et al. 2003 J
: Neurophysiol 89:2707; the same mechanism used in Bazhenov/Timofeev cortical
: slow-oscillation models). No inactivation, no separate time constant --
: I_NaP tracks v instantaneously and is the classic substrate for INTRINSIC
: BURSTING pyramidal cells (L5 "TuftIB", tc_architecture.py): a small
: sustained inward current that, combined with recurrent excitation, turns a
: single spike into a self-sustaining depolarising plateau -- a burst -- that
: SK2/cad (mod/sk2.mod, mod/cad.mod) then terminates. This is what gives a
: cortical column its own robust, periodic slow-oscillation UP-state
: generator, rather than depending on an externally-injected volley for
: every cycle.

NEURON {
    SUFFIX inap
    USEION na READ ena WRITE ina
    RANGE gnabar
}
UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S)  = (siemens)
}
PARAMETER {
    gnabar = 5e-4 (S/cm2)
    vhalf  = -50  (mV)
    slope  = 5    (mV)
}
ASSIGNED {
    v    (mV)
    ena  (mV)
    ina  (mA/cm2)
    minf
}
BREAKPOINT {
    minf = 1 / (1 + exp(-(v - vhalf)/slope))
    ina  = gnabar * minf * (v - ena)
}
