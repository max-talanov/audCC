"""
Auditory input to the thalamocortical model: lemniscal IC -> MGBv (TC) drivers
(PLAN-auditory-input.md, design D1). Stage A: one tonotopic channel, so a
stimulus has a type, onset, duration and level but no frequency.

- IC "fibres" are artificial spike sources, not cells. Each fibre's spike
  train is a deterministic function of (spec seed, fibre id), so every MPI
  rank computes the same trains with no communication, and results do not
  depend on the rank count (the same rule ParallelCorticoThalamicNet follows).
- Each TC cell gets `k` fibres (drawn per TC index), each a strong driver
  synapse with weight `w` (NOT divided by k, unlike the model's k=100
  modulator projections).
- Spikes are delivered with NetCon(None, syn).event(t) from an
  FInitializeHandler, so the same mechanism can later queue closed-loop
  stimuli during a run (design D5).

Rate model per fibre, for a tone of level L starting at t0:

    r(t) = r_spont + A(L) * (1 - exp(-s/rise)) * (sus + (1 - sus) * exp(-s/tau_on))
           for 0 <= s = t - t0 < dur, then an exp(-s'/tau_off) offset decay,
    A(L) = r_max / (1 + exp(-(L - level_50) / level_slope))

drawn as an inhomogeneous Poisson process with an absolute refractory period
(spikes closer than `refractory` to the previous one are dropped). A click
gives each fibre one spike with probability A(L)/r_max, jittered by
`click_jitter` (ms, SD).

The stimulus spec is a JSON file (no YAML dependency, so it runs on MN5 with
only neuron + numpy). See stim/tones_random.json. Either give an explicit
"events" list or a "protocol" that generates one:

    {"protocol": {"type": "random_tones", "start": 2000, "stop": 60000,
                  "isi_min": 2000, "isi_max": 4000, "dur": 50, "level": 60}}
"""

import json

import numpy as np

DEFAULTS = {
    "seed": 1,
    "n_fibres": None,      # None: one fibre per TC cell
    "k": 8,                # fibres per TC cell
    # uS per fibre. aud_checks.py --only tc: a ~4 mV EPSP; at -68 mV two
    # coincident fibres fire the TC cell and one does not; at -80 mV two
    # fibres trigger a 4-spike I_T burst. Calibrated on a depolarised cell
    # only and then held fixed across states (PLAN risk R3).
    "w": 0.002,
    "tau1": 0.5,           # ms, AMPA rise
    "tau2": 3.0,           # ms, AMPA decay (collicular EPSCs decay in a few ms)
    "delay": 1.0,          # ms, IC -> MGB conduction
    "r_spont": 0.0,        # Hz; background input is a Stage B decision
    "r_max": 250.0,        # Hz, peak driven rate
    "level_50": 40.0,      # dB, level at half r_max
    "level_slope": 8.0,    # dB
    "rise": 2.0,           # ms
    "tau_on": 10.0,        # ms, onset-transient decay
    "sus": 0.3,            # sustained / peak ratio
    "tau_off": 3.0,        # ms, offset decay
    "refractory": 1.0,     # ms
    "click_jitter": 0.5,   # ms
    "dt": 0.1,             # ms, rate-model grid
}

TYPE_CODE = {"tone": 0, "click": 1, "noise": 2}


def load_spec(path_or_dict, tstop=np.inf):
    """Spec with defaults filled in; a protocol is generated up to
    min(protocol stop, tstop)."""
    spec = dict(DEFAULTS)
    if isinstance(path_or_dict, dict):
        user = path_or_dict
    else:
        with open(path_or_dict) as f:
            user = json.load(f)
    spec.update({k: v for k, v in user.items() if k not in ("events", "protocol")})
    spec["events"] = list(user.get("events", []))
    if "protocol" in user:
        spec["events"] += make_protocol(user["protocol"], spec["seed"], tstop)
    spec["events"] = sorted((e for e in spec["events"] if e["t"] < tstop),
                            key=lambda e: e["t"])
    return spec


def make_protocol(p, seed, tstop=np.inf):
    """Stimulus schedule. Deterministic in (seed, protocol), independent of
    rank. 'random_tones': uniform random inter-onset intervals, so tones fall
    at all phases of the slow oscillation (Stage C1)."""
    if p["type"] != "random_tones":
        raise ValueError("unknown protocol type %r" % p["type"])
    rng = np.random.default_rng([int(seed), 99])
    kind = p.get("stim", "tone")
    levels = np.atleast_1d(p.get("level", 60.0))
    events, t = [], float(p["start"])
    while t < min(p["stop"], tstop):
        events.append({"t": t, "type": kind, "dur": float(p.get("dur", 50.0)),
                       "level": float(rng.choice(levels))})
        t += rng.uniform(p["isi_min"], p["isi_max"])
    return events


def drive_amplitude(level, spec):
    return spec["r_max"] / (1.0 + np.exp(-(level - spec["level_50"]) / spec["level_slope"]))


def rate_envelope(s, dur, spec):
    """Driven rate shape (0..1 x A) at times s (ms) after onset."""
    s = np.asarray(s, float)
    on = (s >= 0) & (s < dur)
    env = np.zeros_like(s)
    so = s[on]
    env[on] = (1 - np.exp(-so / spec["rise"])) * (
        spec["sus"] + (1 - spec["sus"]) * np.exp(-so / spec["tau_on"]))
    off = s >= dur
    if dur > 0:
        end = (1 - np.exp(-dur / spec["rise"])) * (
            spec["sus"] + (1 - spec["sus"]) * np.exp(-dur / spec["tau_on"]))
        env[off] = end * np.exp(-(s[off] - dur) / spec["tau_off"])
    return env


def fibre_train(fid, spec, tstop):
    """Spike times (ms) of IC fibre `fid`: deterministic in (seed, fid)."""
    rng = np.random.default_rng([int(spec["seed"]), 1, int(fid)])
    dt = spec["dt"]
    spikes = []
    if spec["r_spont"] > 0:
        n = rng.poisson(spec["r_spont"] * tstop / 1000.0)
        spikes.append(rng.uniform(0, tstop, n))
    for ev in spec["events"]:
        A = drive_amplitude(ev["level"], spec)
        if ev["type"] == "click":
            if rng.random() < A / spec["r_max"]:
                spikes.append([ev["t"] + rng.normal(0, spec["click_jitter"])])
            continue
        # tone / noise: one channel in Stage A, so both drive every fibre.
        win = ev["dur"] + 6 * spec["tau_off"]
        s = np.arange(0.0, win, dt)
        p = A * rate_envelope(s, ev["dur"], spec) * dt / 1000.0
        hit = rng.random(len(s)) < p
        spikes.append(ev["t"] + s[hit])
    if not spikes:
        return np.zeros(0)
    t = np.sort(np.concatenate([np.asarray(x, float) for x in spikes]))
    t = t[(t >= 0) & (t < tstop)]
    # absolute refractory period: drop spikes too close to the last kept one
    keep, last = [], -np.inf
    for x in t:
        if x - last >= spec["refractory"]:
            keep.append(x)
            last = x
    return np.asarray(keep)


def fibres_of_tc(tc_index, spec, n_fibres):
    """The k fibre ids converging on TC cell number tc_index (0-based)."""
    rng = np.random.default_rng([int(spec["seed"]), 2, int(tc_index)])
    k = min(spec["k"], n_fibres)
    return rng.choice(n_fibres, size=k, replace=False)


class AuditoryInput:
    """IC -> TC driver input on a ParallelCorticoThalamicNet (or any object
    with .ranges['tc'], .pc, .cells, and ._target_sec)."""

    def __init__(self, net, spec, tstop):
        from neuron import h
        self.h = h
        self.spec = load_spec(spec, tstop)
        self.tstop = tstop
        lo, hi = net.ranges["tc"]
        self.n_fibres = self.spec["n_fibres"] or (hi - lo)
        self.syns, self.ncs = [], []
        trains = {}
        for gid in range(lo, hi):
            if int(net.pc.gid_exists(gid)) == 0:
                continue
            syn = h.Exp2Syn(net.cells[gid].soma(0.5))
            syn.e, syn.tau1, syn.tau2 = 0.0, self.spec["tau1"], self.spec["tau2"]
            self.syns.append(syn)
            for fid in fibres_of_tc(gid - lo, self.spec, self.n_fibres):
                if fid not in trains:
                    trains[fid] = fibre_train(fid, self.spec, tstop)
                nc = h.NetCon(None, syn)
                nc.weight[0] = self.spec["w"]
                self.ncs.append((nc, trains[fid]))
        self._fih = h.FInitializeHandler(self._queue)

    def _queue(self):
        d = self.spec["delay"]
        for nc, train in self.ncs:
            for t in train:
                nc.event(float(t) + d)

    def log(self):
        """Arrays for the run's .npz: the stimulus schedule and every fibre's
        spikes (call on one rank; it recomputes all fibres)."""
        ev = self.spec["events"]
        ids, times = [], []
        for fid in range(self.n_fibres):
            t = fibre_train(fid, self.spec, self.tstop)
            times.append(t)
            ids.append(np.full(len(t), fid))
        return {
            "stim_t": np.array([e["t"] for e in ev], float),
            "stim_type": np.array([TYPE_CODE[e["type"]] for e in ev], int),
            "stim_dur": np.array([e.get("dur", 0.0) for e in ev], float),
            "stim_level": np.array([e["level"] for e in ev], float),
            "ic_times": np.concatenate(times) if times else np.zeros(0),
            "ic_ids": np.concatenate(ids) if ids else np.zeros(0, int),
            "stim_spec": json.dumps({k: v for k, v in self.spec.items() if k != "events"}),
        }
