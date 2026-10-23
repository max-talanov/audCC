"""
Tests for stim_train.py on a synthetic run (cortical activity pulsing at
0.8 Hz = a 1250 ms period; see test_stim_phase._synthetic_npz):

    python3 neuron/test_stim_train.py
"""

import os
import tempfile
import unittest

import numpy as np

import stim_train as ST
import test_stim_phase as TP


def _run_with_schedule(path, sched):
    base = os.path.join(os.path.dirname(path), "base.npz")
    TP._synthetic_npz(base, T=60000.0)
    z = dict(np.load(base, allow_pickle=True))
    z["ol_t"] = np.array([s[0] for s in sched], float)
    z["ol_n"] = np.array([s[1] for s in sched], int)
    z["ol_isi"] = np.array([s[2] for s in sched], float)
    z["ol_dur"] = np.full(len(sched), 50.0)
    z["ol_level"] = np.full(len(sched), 45.0)
    np.savez(path, **z)


class Trains(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        cls.p = os.path.join(cls.d, "run.npz")
        # the train at 1250 ms ISI is in step with the 1250 ms cortical rhythm, the one at 1000 ms drifts
        _run_with_schedule(cls.p, [(10000.0, 5, 1250.0), (30000.0, 5, 1000.0)])

    def test_pulses_of(self):
        z = np.load(self.p, allow_pickle=True)
        P = ST.pulses_of(z)
        self.assertEqual(len(P), 10)
        self.assertEqual([round(x[0]) for x in P[:3]], [10000, 11250, 12500])
        self.assertEqual(P[3][1:], ((5, 1250.0), 0, 3))
        self.assertEqual(P[-1][0], 30000.0 + 4 * 1000.0)

    def test_entrainment_vector_strength_follows_the_rhythm(self):
        pr, tr = ST.run_rows(self.p, skip=2000.0)
        self.assertEqual(len(tr), 2)
        by = {r["isi"]: r for r in tr}
        self.assertGreaterEqual(by[1250.0]["n_up"], 3)
        R1250 = ST.rayleigh(by[1250.0]["phases"])[1]
        R1000 = ST.rayleigh(by[1000.0]["phases"])[1]
        self.assertGreater(R1250, 0.9)
        self.assertLess(R1000, R1250)

    def test_pulse_rows_and_post_rates(self):
        pr, tr = ST.run_rows(self.p, skip=2000.0)
        self.assertEqual({r["k"] for r in pr}, {0, 1, 2, 3, 4})
        self.assertTrue(all(r["n"] == 5 for r in pr))
        for r in tr:
            self.assertGreaterEqual(r["post_up_hz"], 0.0)
            self.assertAlmostEqual(r["span"], 5 * r["isi"])

    def test_summary_and_figures(self):
        pr, tr = ST.run_rows(self.p, skip=2000.0)
        D = {"stim": (pr, tr), "sham": (pr, tr)}
        lines = ST.summary(D)
        self.assertEqual(len(lines), 3)   # header + 2 protocols
        self.assertIn("5 @ 1.00 Hz", lines[1])
        self.assertIn("5 @ 0.80 Hz", lines[2])
        ST.figures(D, os.path.join(self.d, "f"))
        self.assertTrue(os.path.exists(os.path.join(self.d, "f_per_pulse.png")))
        self.assertTrue(os.path.exists(os.path.join(self.d, "f_entrainment.png")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
