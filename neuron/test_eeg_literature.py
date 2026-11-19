import numpy as np
import eeg_literature as L

FS = 1000.0


def test_half_waves_on_sine():
    t = np.arange(0, 60, 1 / FS)
    so = np.sin(2 * np.pi * 0.7 * t)
    p = L.half_waves(so, 1)
    n = L.half_waves(so, -1)
    assert len(p) > 30 and len(n) > 30
    assert abs(p[:, 1].mean() - 1) < 0.01


def test_xcorr_recovers_lag():
    rng = np.random.default_rng(1)
    n = 60000
    a = np.convolve(rng.standard_normal(n), np.ones(300) / 300, "same")
    b = np.roll(a, 400)                                  # b lags a by 0.4 s
    lags, c = L.xcorr(a, b)
    assert abs(lags[np.argmax(c)] - 0.4) < 0.06


def test_spindles_moelle_finds_bursts():
    t = np.arange(0, 30, 1 / FS)
    r = 0.02 + np.zeros_like(t)
    for c in (8.0, 16.0, 24.0):
        r += 0.5 * np.exp(-((t - c) / 0.2) ** 2)
    assert len(L.spindles_moelle(r)) == 3
