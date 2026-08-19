TITLE High-voltage-activated Ca2+ current (cortical pyramidal spike-triggered
: Ca2+ influx)
:
: Reuveni, Friedman, Amitai & Gutnick (1993) J Neurosci 13:4609 HVA kinetics,
: as used in Destexhe/Pospischil-style cortical pyramidal models. No fast
: inactivation -- a broad m^2 activation opening on the upstroke of every
: somatic spike, so Ca2+ influx tracks firing rate. This is the Ca2+ SOURCE for
: the cortical SK2 adaptation current (see sk2.mod / cad.mod): a pyramidal cell
: that fires a burst accumulates submembrane Ca2+, SK2 opens, and the resulting
: after-hyperpolarisation is what ends a cortical UP state -- the same
: SK2 + Ca2+-pool building block already validated in the thalamic RE cell
: (Fernandez & Luthi 2020, sect. V.A.1), reused here for cortical
: spike-frequency adaptation / slow-oscillation bistability (sect. VI).

NEURON {
    SUFFIX ical
    USEION ca READ eca WRITE ica
    RANGE gcabar, m
}
UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S)  = (siemens)
}
PARAMETER {
    gcabar = 1e-4 (S/cm2)
}
ASSIGNED {
    v      (mV)
    eca    (mV)
    ica    (mA/cm2)
    minf
    taum   (ms)
    alpha  (/ms)
    beta   (/ms)
}
STATE { m }
BREAKPOINT {
    SOLVE state METHOD cnexp
    ica = gcabar * m*m * (v - eca)
}
DERIVATIVE state {
    evaluate(v)
    m' = (minf - m)/taum
}
INITIAL {
    evaluate(v)
    m = minf
}
PROCEDURE evaluate(v(mV)) {
    alpha = 0.055 * vtrap(-27 - v, 3.8)
    beta  = 0.94 * exp((-75 - v)/17)
    taum  = 1/(alpha + beta)
    minf  = alpha * taum
}
FUNCTION vtrap(x, y) {
    : x / (exp(x/y) - 1), the standard alpha-rate form (matches hh2.mod's
    : vtrap). NOTE: an earlier version used x/(1-exp(-x/y)) = this * exp(x/y)
    : -- for x/y ~ 11 that is a ~5x10^4 error, which pinned m_inf near 1 at
    : rest (~-68 mV) instead of near 0 and caused spontaneous runaway firing
    : with gcabar as low as 1.2e-4, independent of any bursting mechanism.
    if (fabs(x/y) < 1e-6) { vtrap = y * (1 - x/y/2) }
    else { vtrap = x / (exp(x/y) - 1) }
}
