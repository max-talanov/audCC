"""
Unit tests for stim_phase.py on synthetic data (no MN5 output needed):

    python3 neuron/test_stim_phase.py
"""

import os
import tempfile
import unittest

import numpy as np

import stim_phase as SP

FS = SP.FS


def _osc(f=0.8, T=120.0, noise=0.0, seed=0):
    t = np.arange(0, T, 1.0 / FS)
    x = 3.0 + np.cos(2 * np.pi * f * t)  # cosine: peak at t = 0
    if noise:
        x = x + np.random.default_rng(seed).normal(0, noise, len(t))
    return t, x


class SOPhase(unittest.TestCase):
    def test_bins(self):
        d = np.radians
        self.assertEqual(SP.so_bin(d(-90)), 0)   # rising
        self.assertEqual(SP.so_bin(d(0)), 1)     # peak
        self.assertEqual(SP.so_bin(d(90)), 2)    # falling
        self.assertEqual(SP.so_bin(d(180)), 3)   # trough
        self.assertEqual(SP.so_bin(d(-180)), 3)
        self.assertEqual(SP.so_bin(d(-136)), 3)  # just below the rising bin
        self.assertEqual(SP.so_bin(d(-44)), 1)
        self.assertEqual(SP.so_bin(np.nan), -1)

    def test_acausal_phase_convention(self):
        t, x = _osc()
        ph, amp = SP.phase_acausal(x)
        mid = (t > 20) & (t < 100)
        peaks = np.isclose(((t * 0.8) % 1.0), 0.0, atol=0.5 / (0.8 * FS)) & mid
        troughs = np.isclose(((t * 0.8) % 1.0), 0.5, atol=0.5 / (0.8 * FS)) & mid
        rising = np.isclose(((t * 0.8) % 1.0), 0.75, atol=0.5 / (0.8 * FS)) & mid
        self.assertLess(abs(np.angle(np.mean(np.exp(1j * ph[peaks])))), np.radians(8))
        self.assertGreater(abs(np.angle(np.mean(np.exp(1j * ph[troughs])))), np.radians(172))
        self.assertAlmostEqual(np.degrees(np.angle(np.mean(np.exp(1j * ph[rising])))), -90, delta=8)
        self.assertGreater(amp[mid].mean(), 0.8)

    def test_causal_matches_acausal_on_clean_oscillation(self):
        t, x = _osc(noise=0.1)
        ph, _ = SP.phase_acausal(x)
        errs = []
        for tm in np.arange(20000.0, 100000.0, 733.0):
            pc = SP.phase_causal(x, tm)
            errs.append(abs(np.angle(np.exp(1j * (pc - ph[int(tm // SP.BIN_MS)])))))
        self.assertLess(np.degrees(np.median(errs)), 15)

    def test_causal_uses_only_the_past(self):
        t, x = _osc(noise=0.3)
        tm = 60000.0
        y = x.copy()
        y[int(tm // SP.BIN_MS):] = 50.0 * np.random.default_rng(1).normal(size=len(y) - int(tm // SP.BIN_MS))
        self.assertEqual(SP.phase_causal(x, tm), SP.phase_causal(y, tm))

    def test_causal_needs_history(self):
        _, x = _osc()
        self.assertTrue(np.isnan(SP.phase_causal(x, 3000.0)))


def _volleys(times, part=0.5):
    return [{"t": float(t), "re_part": part + 0.01 * k, "tc_part": 0.3, "re_sd": 3.0}
            for k, t in enumerate(times)]


class Spindles(unittest.TestCase):
    def test_trains_need_min_cycles(self):
        v = _volleys([1000 + 80 * k for k in range(8)] + [5000, 5080, 5160])
        tr = SP.trains_of(v, min_cycles=6)
        self.assertEqual(len(tr), 1)
        self.assertEqual(tr[0]["cycles"], 8)
        self.assertEqual(tr[0]["t_peak"], 1000 + 80 * 7)  # re_part grows with k

    def test_gap_splits_trains(self):
        v = _volleys([1000 + 80 * k for k in range(6)] + [2000 + 80 * k for k in range(6)])
        self.assertEqual(len(SP.trains_of(v, 6)), 2)

    def test_spindle_state_is_causal(self):
        v = _volleys([1000 + 80 * k for k in range(8)])
        self.assertFalse(SP.spindle_state(v, 900)[0])
        self.assertFalse(SP.spindle_state(v, 1100)[0])      # only 2 volleys so far
        pres, tp = SP.spindle_state(v, 1000 + 80 * 4 + 10)  # 5 volleys so far
        self.assertTrue(pres)
        self.assertEqual(tp, 1000 + 80 * 4)                 # peak *so far*, not the final one
        self.assertTrue(SP.spindle_state(v, 1000 + 80 * 7 + 100)[0])
        self.assertFalse(SP.spindle_state(v, 1000 + 80 * 7 + 400)[0])


def _synthetic_npz(path, T=40000.0, seed=0):
    """Cortical populations modulated at 0.8 Hz (peak at t = 0 mod 1250 ms), RE
    volley train of 8 volleys at 12.5 Hz starting at 10.0 s."""
    rng = np.random.default_rng(seed)
    R = {"tc": (0, 20), "re": (20, 60), "l4e": (60, 100), "l4i": (100, 110),
         "l23e": (110, 150), "l23i": (150, 160), "l5e": (160, 200), "l5i": (200, 210),
         "l6e": (210, 250), "l6i": (250, 260)}
    ts, gs = [], []
    grid = np.arange(0.0, T, 1.0)
    rate = 1.0 + 40.0 * (0.5 * (1 + np.cos(2 * np.pi * 0.8 * grid / 1000.0))) ** 4  # Hz/cell
    for name in ("l4e", "l23e", "l5e", "l6e"):
        lo, hi = R[name]
        for gid in range(lo, hi):
            n = rng.random(len(grid)) < rate / 1000.0
            ts.append(grid[n])
            gs.append(np.full(n.sum(), gid))
    for k in range(8):
        for gid in range(20, 60):
            ts.append(np.array([10000.0 + 80 * k + rng.uniform(-3, 3)]))
            gs.append(np.array([gid]))
    t, g = np.concatenate(ts), np.concatenate(gs)
    o = np.argsort(t)
    np.savez(path, times=t[o], gids=g[o], ranges=np.array(R, dtype=object), tstop=T,
             sizes=np.array(0))


class RunClass(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        cls.p = os.path.join(cls.d, "syn.npz")
        _synthetic_npz(cls.p)
        cls.sr = SP.Run(cls.p)

    def test_phase_at_peaks(self):
        # activity peaks at multiples of 1250 ms
        ph = [self.sr.phase_a[self.sr.at(1250.0 * k)] for k in range(8, 28)]
        self.assertLess(abs(np.degrees(np.angle(np.mean(np.exp(1j * np.array(ph)))))), 25)

    def test_classify_bins_and_causality(self):
        r = self.sr.classify(20000.0 - 1250 * 0.25)   # a quarter cycle before a peak = rising
        self.assertEqual(r["so_bin"], 0)
        r = self.sr.classify(20000.0 + 1250 * 0.5)    # trough
        self.assertEqual(r["so_bin"], 3)

    def test_spindle_present_and_phase(self):
        r = self.sr.classify(10000.0 + 80 * 5 + 10)
        self.assertTrue(r["sp_present"])
        self.assertGreater(r["sp_bin"], 0)
        self.assertFalse(self.sr.classify(25000.0)["sp_present"])
        self.assertEqual(self.sr.classify(25000.0)["sp_bin"], 0)

    def test_respond(self):
        hit = self.sr.respond(9800.0)     # spindle peak (volley 8 region) inside 1 s
        miss = self.sr.respond(25000.0)
        self.assertTrue(hit["sp_any"])
        self.assertFalse(miss["sp_any"])
        self.assertGreater(hit["sp_pow"], miss["sp_pow"] * 5)
        self.assertGreater(hit["so_amp"], 0)

    def test_tone_rows_skip_edges(self):
        rows = SP.tone_rows(self.sr, [500.0, 15000.0, 39500.0], skip=2000.0)
        self.assertEqual([r["t"] for r in rows], [15000.0])


class Stats(unittest.TestCase):
    @staticmethod
    def _rows(mean, n, seed_ids=(0, 1, 2), so=1, sp=0, sd=1.0, rng=None):
        rng = rng or np.random.default_rng(0)
        return [{"so_bin": so, "sp_bin": sp, "seed": k, "x": float(rng.normal(mean, sd))}
                for k in seed_ids for _ in range(n)]

    def test_true_effect_inside_ci_and_null_includes_zero(self):
        rng = np.random.default_rng(3)
        stim = self._rows(1.0, 60, rng=rng)
        sham = self._rows(0.0, 60, rng=rng)
        e, lo, hi, na, nb = SP.bootstrap_effect(stim, sham, "x", 1, 0, 1000)
        self.assertEqual((na, nb), (180, 180))
        self.assertLess(lo, 1.0)
        self.assertGreater(hi, 1.0)
        self.assertGreater(lo, 0.0)
        null = self._rows(0.0, 60, rng=rng)
        e, lo, hi, *_ = SP.bootstrap_effect(null, sham, "x", 1, 0, 1000)
        self.assertLess(lo, 0.0)
        self.assertGreater(hi, 0.0)

    def test_bootstrap_is_deterministic_and_nan_for_empty_cell(self):
        rng = np.random.default_rng(4)
        a, b = self._rows(1.0, 20, rng=rng), self._rows(0.0, 20, rng=rng)
        self.assertEqual(SP.bootstrap_effect(a, b, "x", 1, 0, 300, seed=7),
                         SP.bootstrap_effect(a, b, "x", 1, 0, 300, seed=7))
        self.assertTrue(np.isnan(SP.bootstrap_effect(a, b, "x", 2, 0, 100)[0]))

    def test_more_seeds_narrow_ci(self):
        rng = np.random.default_rng(5)
        w = lambda ids: SP.bootstrap_effect(self._rows(1, 30, ids, rng=rng), self._rows(0, 30, ids, rng=rng),
                                            "x", 1, 0, 800)
        few, many = w((0, 1, 2)), w(tuple(range(12)))
        self.assertGreater(few[2] - few[1], many[2] - many[1])

    def test_circular_summary(self):
        rng = np.random.default_rng(6)
        m, R, p, n = SP.circ_summary(rng.normal(0.7, 0.3, 200))
        self.assertAlmostEqual(m, 0.7, delta=0.1)
        self.assertGreater(R, 0.9)
        self.assertLess(p, 1e-6)
        m, R, p, n = SP.circ_summary(rng.uniform(-np.pi, np.pi, 200))
        self.assertGreater(p, 0.01)
        self.assertEqual(SP.circ_summary([np.nan, np.nan])[3], 0)
        self.assertEqual(SP.circ_summary([0.1, np.nan, 0.2])[3], 2)

    def test_grid_counts(self):
        rows = [{"so_bin": 0, "sp_bin": 0, "x": 1.0}, {"so_bin": 0, "sp_bin": 0, "x": 3.0},
                {"so_bin": 2, "sp_bin": 3, "x": 5.0}]
        G, N = SP.grid(rows, "x")
        self.assertEqual(N[0, 0], 2)
        self.assertEqual(G[0, 0], 2.0)
        self.assertEqual(N[2, 3], 1)
        self.assertTrue(np.isnan(G[1, 1]))
        self.assertEqual(N.sum(), 3)


class Schedule(unittest.TestCase):
    kw = dict(gap=(5000.0, 9000.0), n_choices=(1, 2), isi=1075.0, durs=(25.0, 50.0, 100.0, 200.0),
              levels=(30.0, 45.0, 60.0, 75.0))

    def test_deterministic_and_seed_dependent(self):
        import closed_loop as CL
        a = CL.tone_schedule(1, 400000.0, **self.kw)
        self.assertEqual(a, CL.tone_schedule(1, 400000.0, **self.kw))
        self.assertNotEqual(a, CL.tone_schedule(2, 400000.0, **self.kw))

    def test_spacing_bounds_and_choices(self):
        import closed_loop as CL
        T = 400000.0
        sch = CL.tone_schedule(3, T, **self.kw)
        self.assertGreater(len(sch), 40)
        self.assertGreaterEqual(sch[0]["t"], 2000.0 + 5000.0)
        for a, b in zip(sch, sch[1:]):
            end = a["t"] + (a["n"] - 1) * a["isi"]
            self.assertGreaterEqual(b["t"] - end, 5000.0 - 1e-6)
            self.assertLessEqual(b["t"] - end, 9000.0 + 1e-6)
        last = sch[-1]
        self.assertLessEqual(last["t"] + (last["n"] - 1) * last["isi"] + 1500.0, T)
        for x in sch:
            self.assertIn(x["n"], (1, 2))
            self.assertIn(x["dur"], self.kw["durs"])
            self.assertIn(x["level"], self.kw["levels"])
        self.assertEqual({x["n"] for x in sch}, {1, 2})

    def test_factors_are_balanced_over_many_tones(self):
        import closed_loop as CL
        lv = [x["level"] for sd in range(12) for x in CL.tone_schedule(sd, 400000.0, **self.kw)]
        for v in self.kw["levels"]:
            self.assertGreater(lv.count(v) / len(lv), 0.18)

    def test_isi_list_and_scalar_compat(self):
        import closed_loop as CL
        kw = dict(self.kw)
        kw.pop("isi")
        a = CL.tone_schedule(7, 400000.0, isi=1075.0, **kw)
        b = CL.tone_schedule(7, 400000.0, isi=[1075.0], **kw)
        self.assertEqual(a, b)
        many = CL.tone_schedule(7, 400000.0, isi=[600.0, 1000.0, 1250.0, 2000.0],
                                gap=(5000.0, 9000.0), n_choices=(2, 5, 10), durs=(50.0,), levels=(35.0,))
        self.assertEqual({x["isi"] for x in many}, {600.0, 1000.0, 1250.0, 2000.0})
        for x, y in zip(many, many[1:]):
            gap = y["t"] - (x["t"] + (x["n"] - 1) * x["isi"])
            self.assertGreaterEqual(gap, 5000.0 - 1e-6)
            self.assertLessEqual(gap, 9000.0 + 1e-6)

    def test_shorter_run_is_a_prefix(self):
        import closed_loop as CL
        long = CL.tone_schedule(5, 400000.0, **self.kw)
        short = CL.tone_schedule(5, 100000.0, **self.kw)
        self.assertEqual(short, long[:len(short)])


class ToneInfo(unittest.TestCase):
    def test_open_loop_and_closed_loop_tone_sources(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "ol.npz")
        np.savez(p, ol_t=np.array([8000.0, 20000.0]), ol_n=np.array([1, 2]),
                 ol_dur=np.array([50.0, 100.0]), ol_level=np.array([45.0, 60.0]),
                 ol_isi=np.array([1075.0, 1075.0]), tstop=40000.0)
        t, info = SP.tones_of_run(np.load(p, allow_pickle=True))
        self.assertEqual(list(t), [8000.0, 20000.0])
        self.assertEqual(info[1], {"n": 2, "dur": 100.0, "level": 60.0, "isi": 1075.0})
        p = os.path.join(d, "cl.npz")
        np.savez(p, cl_t_det=np.array([3000.0, 6000.0]), cl_d=np.array([200.0, 400.0]),
                 cl_stim=np.array([1, 0]), tstop=40000.0)
        t, info = SP.tones_of_run(np.load(p, allow_pickle=True))
        self.assertEqual(list(t), [3200.0])
        self.assertEqual(info, [{"d": 200.0}])


if __name__ == "__main__":
    unittest.main(verbosity=2)
