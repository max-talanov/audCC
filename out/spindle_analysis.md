# Spindle analysis (per layer)

Auditory thalamo-cortical column, 15 s, neuron model `iaf_cond_exp`. Sigma band 10-15 Hz.

| layer | region | freq (Hz) | density /min | dur (s) | rel. power | lag vs MGB (ms) | SO-coupled |
|---|---|---|---|---|---|---|---|
| MGB | thalamus | 13.0 | 60.0 | 0.46 | 1.00 | 0 | 0% |
| nRT | thalamus | 13.0 | 60.0 | 0.44 | 4.20 | 0 | 0% |
| L4 | cortex | 13.0 | 60.0 | 0.45 | 5.15 | 11 | 0% |
| L23 | cortex | 13.0 | 120.0 | 0.22 | 3.02 | 68 | 50% |
| L5 | cortex | 13.0 | 120.0 | 0.22 | 3.68 | 85 | 50% |
| L6 | cortex | 13.0 | 120.0 | 0.22 | 1.24 | 94 | 50% |

**Thalamus -> cortex propagation:** the sigma-band spindle envelope in the cortical layers lags MGB by 64 ms on average (positive = spindle arrives in cortex after the thalamus) -- the spindle travelling up the column.

> Notes: frequency, relative power and the propagation lag are the robust measures. Density reflects the imposed-drive `iaf_cond_exp` regime (~one spindle per 1 Hz slow-oscillation cycle, i.e. ~60/min), higher than the 2-8/min of natural NREM; the per-layer detection threshold also splits some cortical spindles, inflating the upper-layer count. Use the AdEx/HH configs for physiological spindle *density* (see the validator, `tc_validate.py`).
