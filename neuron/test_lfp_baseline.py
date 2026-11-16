import unittest

import numpy as np

import lfp_baseline as LB


class TestSpectra(unittest.TestCase):
    def test_band_fractions_sum_to_one(self):
        fs = 1000.0
        t = np.arange(0, 60, 1 / fs)
        x = np.sin(2 * np.pi * 0.8 * t) + 0.5 * np.sin(2 * np.pi * 13 * t)
        f, p = LB.psd(x, fs)
        fr = LB.band_fractions(f, p)
        self.assertAlmostEqual(sum(fr.values()), 1.0, delta=0.05)
        self.assertGreater(fr["SO 0.3-1.5"], 0.7)
        self.assertGreater(fr["sigma 9-16"], 0.05)

    def test_peak_freq(self):
        fs = 200.0
        t = np.arange(0, 60, 1 / fs)
        f, p = LB.psd(np.sin(2 * np.pi * 12.0 * t), fs)
        self.assertAlmostEqual(LB.peak_freq(f, p, 5.0, 30.0), 12.0, delta=0.2)

    def test_lfp_traces_names_and_flip(self):
        rng = np.random.default_rng(0)
        R = {"tc": (0, 10), "re": (10, 20), "l4e": (20, 30), "l4i": (30, 40), "l23e": (40, 50),
             "l23i": (50, 60), "l5e": (60, 70), "l5i": (70, 80), "l6e": (80, 90), "l6i": (90, 100)}
        T = 20000.0
        t = np.sort(rng.uniform(0, T, 5000))
        g = rng.integers(0, 100, 5000)
        z = {"times": t, "gids": g, "ranges": np.array(R, dtype=object), "tstop": np.array(T)}

        class Z(dict):
            def __getitem__(self, k):
                v = dict.__getitem__(self, k)
                return v
        z = Z(times=t, gids=g, ranges=np.array(R, dtype=object), tstop=np.array(T))
        L = LB.lfp_traces(z)
        for k in LB.GROUPS + ["L2/3", "L4", "L5", "L6", "cortex", "thalamus", "network"]:
            self.assertIn(k, L)
            self.assertEqual(len(L[k]), int(T))
        self.assertAlmostEqual(float(np.std(L["tc"])), 1.0, delta=0.01)


if __name__ == "__main__":
    unittest.main()
