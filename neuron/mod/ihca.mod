: Hyperpolarisation-activated cation current I_h with CALCIUM-DEPENDENT
: up-regulation -- Destexhe, Bal, McCormick & Sejnowski (1996), J Neurophysiol
: 76:2049; Destexhe et al. (1993) Biophys J 65:1538.
:
: This is the canonical Ca2+-dependent thalamic spindle mechanism. I_h opens on
: hyperpolarisation (C <-> O), but the open form is additionally stabilised by a
: Ca2+-bound regulating factor:
:
:     C   <-> O                    voltage-dependent opening (alpha, beta)
:     P0 + nca Ca  <-> P1          Ca2+ binds the regulating factor (k1ca, k2)
:     O  + P1      <-> OL          the "locked open" form         (k3p, k4)
:
:     I_h = ghbar * (O + ginc*OL) * (v - eh)
:
: Ca2+ entering through I_T on successive rebound bursts accumulates, shifts the
: equilibrium toward the locked-open OL state and so progressively up-regulates
: I_h. The resulting depolarisation removes I_T de-inactivation and TERMINATES
: the spindle; the slow unbinding (k2 = 4e-4 /ms, tau ~ 2.5 s) then produces the
: 5-10 s refractory period before the next spindle (Fernandez & Luthi 2020,
: sect. V.D.1). Without this loop a thalamic network expresses the slow intrinsic
: delta-like oscillation rather than waxing/waning spindles.
:
: Ca2+ SOURCE. This reads the submembrane pool maintained by cad.mod, which lives
: in the private "sk" ion species -- so the Ca2+ that regulates I_h never feeds
: back into the T-current's GHK reversal potential. Insert cad alongside this
: mechanism or the pool stays at rest and I_h reduces to plain voltage-gated I_h.
:
: The states are nonlinearly coupled (k1ca depends on Ca, k3p on p1), so this
: integrates with derivimplicit, NOT cnexp.

NEURON {
    SUFFIX ihca
    NONSPECIFIC_CURRENT i
    USEION sk READ ski VALENCE 2
    RANGE ghbar, eh, shift, h_inf, tau_s, ginc, cac, k2, k4, Pc, i
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S)  = (siemens)
}

PARAMETER {
    ghbar = 2e-5   (S/cm2)
    eh    = -43    (mV)
    shift = 0      (mV)
    ginc  = 2              : the locked-open form conducts ginc-fold
    cac   = 0.002          : half-activation Ca2+ of the regulating factor
    Pc    = 0.01           : half-activation of the O + P1 -> OL step
    k2    = 0.0004 (/ms)   : Ca2+ unbinding  -- sets the ~2.5 s refractory decay
    k4    = 0.001  (/ms)   : OL -> O + P1
    nca   = 4              : Ca2+ ions bound per regulating factor
    nexp  = 1
    taum  = 20     (ms)
    celsius (degC)
}

ASSIGNED {
    v      (mV)
    ski
    i      (mA/cm2)
    h_inf
    tau_s  (ms)
    alpha  (/ms)
    beta   (/ms)
    k1ca   (/ms)
    k3p    (/ms)
    o1
    p0
}

STATE { c1 o2 p1 }

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    o1 = 1 - c1 - o2
    i  = ghbar * (o1 + ginc*o2) * (v - eh)
}

DERIVATIVE states {
    evaluate(v, ski)
    o1 = 1 - c1 - o2
    p0 = 1 - p1
    c1' = beta*o1  - alpha*c1
    p1' = k1ca*p0  - k2*p1
    o2' = k3p*o1   - k4*o2
}

INITIAL {
    evaluate(v, ski)
    : steady state with negligible Ca2+ binding
    c1 = 1 - h_inf
    o2 = 0
    p1 = 0
}

PROCEDURE evaluate(v(mV), ski) {
    LOCAL vs, tadj
    vs   = v + shift
    tadj = 3 ^ ((celsius - 36)/10)
    h_inf = 1 / (1 + exp((vs + 75)/5.5))
    tau_s = (taum + 1000 / (exp((vs + 71.5)/14.2) + exp(-(vs + 89)/11.6))) / tadj
    alpha = h_inf / tau_s
    beta  = (1 - h_inf) / tau_s
    : Ca2+-dependent forward rates (Destexhe et al. 1993, eqs. 6-7)
    k1ca = k2 * (ski/cac) ^ nca
    k3p  = k4 * (p1 /Pc ) ^ nexp
}
