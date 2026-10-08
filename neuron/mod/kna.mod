TITLE Na+-activated K+ current (K_Na, Slack/Slick) with a slow Na+ pool

COMMENT
TRN-side spindle termination and refractoriness (Fernandez & Luthi 2020,
sect. V.D.1 and Fig. 7g): during repeated burst cycles, Na+ entering with
the spikes accumulates and activates K_Na channels. The TRN cell then
hyperpolarises and its conductance rises, so later bursts have fewer spikes
and TC EPSPs are shunted. After the spindle, the K_Na current decays with
the Na+ pool over seconds, which is the TRN refractory period.

Phenomenological, after the K_Na of Compte et al. 2003 / Bazhenov et al.
2002, in a form that is easy to calibrate:

    nax' = alpha * (-ina) - nax / taur        excess Na+ in a lumped pool
    w    = nax^nh / (nax^nh + kd^nh)          K_Na activation, 0..1
    ik   = gbar * w * (v - ek)

nax is the EXCESS Na+ over rest, in a lumped submembrane pool. alpha is
calibrated, not derived from a physical shell: it is set so that one spindle
of RE bursts raises nax to about kd (w ~ 0.5); see PLAN-spindels.md. The pool
is private (STATE of this mechanism), so ena and the spikes are unchanged.
It READs ina from hh2 (the Na+ source), like cad.mod reads ica.
ENDCOMMENT

NEURON {
    SUFFIX kna
    USEION na READ ina
    USEION k READ ek WRITE ik
    RANGE gbar, alpha, taur, kd, nh, nax, w, ik
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mM) = (milli/liter)
}

PARAMETER {
    gbar  = 0      (S/cm2)
    alpha = 2      (mM cm2/mA ms)     : lumped gain: an RECell spike (0.25 mA ms/cm2 of Na+) adds 0.5 mM, ~20 spikes per spindle reach kd
    taur  = 4000   (ms)               : Na+ removal (pump), sets refractoriness
    kd    = 10     (mM)               : half-activation of K_Na by excess Na+
    nh    = 3.5                       : Hill coefficient (Compte et al. 2003)
}

ASSIGNED {
    v   (mV)
    ek  (mV)
    ina (mA/cm2)
    ik  (mA/cm2)
    w
}

STATE { nax (mM) }

INITIAL {
    nax = 0
    w = 0
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    w = act(nax)
    ik = gbar * w * (v - ek)
}

DERIVATIVE state {
    LOCAL drive
    drive = -alpha * ina
    if (drive < 0) { drive = 0 }
    nax' = drive - nax / taur
}

FUNCTION act(x (mM)) {
    if (x <= 0) {
        act = 0
    } else {
        act = 1 / (1 + pow(kd / x, nh))
    }
}
