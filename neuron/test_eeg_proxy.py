import numpy as np
import eeg_proxy as P

FS = 1000.0


def test_so_detection_on_sine():
    t = np.arange(0, 60, 1 / FS)
    x = np.sin(2 * np.pi * 0.7 * t)
    idx, amp, dur, *_ = P.detect_so(x)
    assert len(idx) > 5 and abs(dur.mean() - 1 / 0.7) < 0.1


def test_spindle_detection_on_burst():
    rng = np.random.default_rng(0)
    t = np.arange(0, 30, 1 / FS)
    x = 0.05 * rng.standard_normal(len(t))
    for c in (8.0, 16.0, 24.0):
        m = np.abs(t - c) < 0.4
        x[m] += np.sin(2 * np.pi * 13 * t[m]) * np.hanning(m.sum())
    so = P.bandpass(x, 0.16, 1.25, 2)
    sp, *_ = P.detect_spindles(x, so)
    assert len(sp) == 3 and abs(sp[:, 3].mean() - 13) < 1


def test_head_filter_removes_dc():
    n = 30000
    x = np.ones(n) + np.sin(2 * np.pi * 2 * np.arange(n) / FS)
    assert abs(P.head_filter(x)[n // 4:3 * n // 4].mean()) < 0.05
