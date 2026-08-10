"""
HDF5 I/O for simulation results.

Lets a cluster run (MareNostrum 5) save its *raw* output once, so every figure,
validation and analysis can be regenerated locally without re-simulating:

    MN5:    tc_run.py ... --save-h5 out/run.h5      # simulate, save raw data
    local:  scp ...:out/run.h5 .
            tc_plot_h5.py run.h5 --outdir out       # all figures + validation

The file stores exactly the ``tc_validate`` result contract (spikes / traces /
meta), so a loaded file can be passed straight to ``tc_validate.validate``,
``tc_run.make_plot``, ``tc_analyze``, etc.

Layout::

    /                       attrs: tstop, dt, seed, neuron_model, config,
                                   n_per_layer (JSON), created, format_version
    /spikes/<layer>/times     float64[]   spike times, ms
    /spikes/<layer>/senders   int64[]     per-spike cell id
    /traces/<layer>/time      float64[]   sample times, ms
    /traces/<layer>/voltage   float64[]   layer-mean V_m, mV

Datasets are gzip-compressed; a 200 s 5k-neuron run is typically a few tens of MB.
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

FORMAT_VERSION = 1


def save_h5(spikes, traces, meta, path, compression="gzip", compression_opts=4):
    """Write a simulation result to HDF5. Returns the path."""
    import h5py

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.attrs["format_version"] = FORMAT_VERSION
        f.attrs["created"] = datetime.now(timezone.utc).isoformat()
        for k, v in (meta or {}).items():
            if isinstance(v, dict):
                f.attrs[k] = json.dumps(v)          # e.g. n_per_layer
            elif v is None:
                continue
            else:
                f.attrs[k] = v

        gs = f.create_group("spikes")
        for layer, d in (spikes or {}).items():
            g = gs.create_group(layer)
            g.create_dataset("times", data=np.asarray(d.get("times", []), float),
                             compression=compression,
                             compression_opts=compression_opts)
            g.create_dataset("senders",
                             data=np.asarray(d.get("senders", []), np.int64),
                             compression=compression,
                             compression_opts=compression_opts)

        gt = f.create_group("traces")
        for layer, d in (traces or {}).items():
            g = gt.create_group(layer)
            g.create_dataset("time", data=np.asarray(d.get("time", []), float),
                             compression=compression,
                             compression_opts=compression_opts)
            g.create_dataset("voltage",
                             data=np.asarray(d.get("voltage", []), float),
                             compression=compression,
                             compression_opts=compression_opts)
    return path


def load_h5(path):
    """Read a result written by :func:`save_h5`. Returns (spikes, traces, meta)."""
    import h5py

    spikes, traces, meta = {}, {}, {}
    with h5py.File(Path(path), "r") as f:
        for k, v in f.attrs.items():
            if isinstance(v, (bytes, np.bytes_)):
                v = v.decode()
            if k == "n_per_layer" and isinstance(v, str):
                try:
                    v = json.loads(v)
                except json.JSONDecodeError:
                    pass
            meta[k] = v
        if "spikes" in f:
            for layer in f["spikes"]:
                g = f["spikes"][layer]
                spikes[layer] = {"times": g["times"][:], "senders": g["senders"][:]}
        if "traces" in f:
            for layer in f["traces"]:
                g = f["traces"][layer]
                traces[layer] = {"time": g["time"][:], "voltage": g["voltage"][:]}
    # tstop is required downstream; recover it if an older file lacks it
    if "tstop" not in meta:
        ends = [t["time"][-1] for t in traces.values() if len(t["time"])]
        ends += [s["times"].max() for s in spikes.values() if len(s["times"])]
        meta["tstop"] = float(max(ends)) if ends else 0.0
    meta.setdefault("seed", 0)
    return spikes, traces, meta


def describe(path):
    """Human-readable summary of an HDF5 result file."""
    spikes, traces, meta = load_h5(path)
    size_mb = Path(path).stat().st_size / 1e6
    lines = [f"{path}  ({size_mb:.1f} MB)",
             f"  tstop={meta.get('tstop', 0)/1000:.1f} s  "
             f"model={meta.get('neuron_model', '?')}  seed={meta.get('seed', '?')}"]
    cfg = meta.get("config")
    if cfg:
        lines.append(f"  config={cfg}")
    npl = meta.get("n_per_layer", {})
    lines.append("  spikes:")
    for layer, d in spikes.items():
        n = len(d["times"])
        cells = npl.get(layer) if isinstance(npl, dict) else None
        rate = (1000.0 * n / (meta["tstop"] * cells)) if cells and meta.get("tstop") else None
        lines.append(f"    {layer:<5} {n:>9,} spikes"
                     + (f"  ({cells} cells, {rate:5.1f} Hz/neuron)" if rate else ""))
    if traces:
        lines.append(f"  traces: {', '.join(sorted(traces))}")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        raise SystemExit("usage: python3 tc_io.py <result.h5>")
    for p in sys.argv[1:]:
        print(describe(p))
