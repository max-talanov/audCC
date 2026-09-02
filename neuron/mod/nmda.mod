TITLE NMDA synapse with Mg2+ block (Jahr & Stevens 1990)

COMMENT
Double-exponential conductance (rise tau1, decay tau2) gated by the
classic Mg2+-block voltage dependence (Jahr & Stevens 1990). Unlike
Exp2Syn (used everywhere else in this model), the net conductance here
is NONLINEAR in voltage: weak near rest (Mg2+ blocks the channel),
strongly unblocking as the cell depolarises. That saturating
nonlinearity is what a stable, bounded recurrent-excitation UP state
needs and a linear synapse cannot provide -- see ctx_thalamus_mpi.py's
_wire_l5_recurrent Option A/B comparison (the plain-Exp2Syn version,
mechanism="exp2syn", had essentially no stable middle ground between
"no effect" and "runaway").
ENDCOMMENT

NEURON {
    POINT_PROCESS NMDA
    RANGE tau1, tau2, e, i, mg, g
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (mM) = (milli/liter)
}

PARAMETER {
    tau1 = 5   (ms)  : rise
    tau2 = 100 (ms)  : decay -- slow, unlike every other synapse (2-8ms)
    e = 0      (mV)
    mg = 1     (mM)  : extracellular Mg2+ (1 mM is physiological)
}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
    factor
}

STATE {
    A (uS)
    B (uS)
}

INITIAL {
    LOCAL tp
    if (tau1/tau2 > 0.9999) {
        tau1 = 0.9999*tau2
    }
    A = 0
    B = 0
    tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
    factor = -exp(-tp/tau1) + exp(-tp/tau2)
    factor = 1/factor
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    g = (B - A) * mgblock(v)
    i = g*(v - e)
}

DERIVATIVE state {
    A' = -A/tau1
    B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
    A = A + weight*factor
    B = B + weight*factor
}

FUNCTION mgblock(v (mV)) {
    : Jahr & Stevens 1990: 1 / (1 + [Mg2+]*exp(-0.062*v)/3.57)
    mgblock = 1 / (1 + exp(-0.062(/mV)*v) * (mg/3.57(mM)))
}
