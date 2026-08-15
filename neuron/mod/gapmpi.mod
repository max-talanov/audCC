: Gap junction for MPI (cross-rank) electrical coupling.
:
: Identical physics to gap.mod, with ONE difference: vgap is a RANGE variable
: rather than a POINTER. ParallelContext.target_var() writes the remote cell's
: membrane potential into vgap every dt, and it can only write into a real
: variable -- taking _ref_ of an unassigned POINTER raises "Invalid pointer".
:
: gap.mod (POINTER + h.setpointer) stays the serial path; this is the parallel
: one. Same g, same current, so the two give the same coupling.

NEURON {
    POINT_PROCESS GapMPI
    RANGE g, i, vgap
    NONSPECIFIC_CURRENT i
}

PARAMETER { g = 0 (uS) }

ASSIGNED {
    v    (mV)
    vgap (mV)
    i    (nA)
}

INITIAL { vgap = v }

BREAKPOINT { i = g * (v - vgap) }
