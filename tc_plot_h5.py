"""
Regenerate every figure, validation and analysis from a saved HDF5 result.

Intended workflow: simulate once on the cluster, plot as often as you like
locally -- no NEST needed here, and no re-simulation.

    # on MareNostrum 5
    sbatch --export=ALL,CONFIG=config/network_auditory_mn5_5k.yaml,TSTOP=200000,SAVE_H5=1 run.sh

    # locally
    scp USER@glogin1.bsc.es:'~/audCC/out/*.h5' .
    python3 tc_plot_h5.py run.h5 --outdir out

Produces (tagged with the file stem):
    <tag>.png            5-panel overview (raster, LFP, zoom, spectrogram, PSD)
    <tag>_decomp.png     slow-wave + spindle band decomposition
    <tag>_layers.png     the spindle propagating up the auditory column
    <tag>_present.png    presentation figure (slow waves + spindles)
    printed             Fernandez & Luthi (2020) validation table + per-layer metrics

Requires only numpy / scipy / matplotlib / h5py.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tc_io                      # noqa: E402
import tc_run                     # noqa: E402


def _rhythms(spikes, tstop):
    """Slow-wave and spindle peaks, as tc_run.main computes them."""
    cort = [l for l in ["L23", "L5", "L6"] if l in spikes]
    thal = [l for l in ["MGB", "nRT"] if l in spikes]
    slow_peak = spindle_peak = 0.0
    if cort:
        allc = np.concatenate([spikes[l]["times"] for l in cort])
        _, rc = tc_run.population_rate(allc, tstop, bin_ms=5.0, smooth_ms=20.0)
        slow_peak, _, _, _ = tc_run.detect_peak(rc[100:], 200.0, 0.3, 1.5)
    if thal:
        allt = np.concatenate([spikes[l]["times"] for l in thal])
        _, rt = tc_run.population_rate(allt, tstop, bin_ms=1.0, smooth_ms=3.0)
        spindle_peak, _, _, _ = tc_run.detect_peak(rt[500:], 1000.0, 9.0, 16.0)
    return slow_peak, spindle_peak


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5", help="HDF5 result written by tc_run.py --save-h5")
    ap.add_argument("--outdir", type=str, default="out")
    ap.add_argument("--tag", type=str, default=None,
                    help="output filename tag (default: the .h5 file stem)")
    ap.add_argument("--no-validate", action="store_true",
                    help="skip the Fernandez & Luthi criteria table")
    ap.add_argument("--no-analyze", action="store_true",
                    help="skip the per-layer metrics table")
    ap.add_argument("--info", action="store_true",
                    help="print a summary of the file and exit")
    args = ap.parse_args(argv)

    if args.info:
        print(tc_io.describe(args.h5))
        return 0

    print(tc_io.describe(args.h5))
    spikes, traces, meta = tc_io.load_h5(args.h5)
    tstop = meta["tstop"]
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    tag = args.tag or Path(args.h5).stem

    slow_peak, spindle_peak = _rhythms(spikes, tstop)
    print(f"\n  slow-wave peak : {slow_peak:.2f} Hz")
    print(f"  spindle peak   : {spindle_peak:.2f} Hz")

    # ---- figures ----
    print("\nfigures:")
    try:
        tc_run.make_plot(spikes, traces, meta, slow_peak, spindle_peak,
                         outdir / f"{tag}.png")
    except Exception as e:
        print(f"  (overview failed: {e})")
    try:
        tc_run.make_decomposition_plot(traces, meta, slow_peak, spindle_peak,
                                       outdir / f"{tag}_decomp.png")
    except Exception as e:
        print(f"  (decomposition failed: {e})")
    try:
        tc_run.make_layer_spindle_plot(spikes, meta, outdir / f"{tag}_layers.png")
    except Exception as e:
        print(f"  (per-layer failed: {e})")
    try:
        import tc_present
        tc_present.make_figure(spikes, traces, meta,
                               outdir / f"{tag}_present.png")
    except Exception as e:
        print(f"  (presentation figure failed: {e})")

    # ---- validation ----
    if not args.no_validate:
        try:
            import tc_validate
            tc_validate.validate_result(spikes, traces, meta)
        except Exception as e:
            print(f"(validation failed: {e})")

    # ---- per-layer analysis ----
    if not args.no_analyze:
        try:
            import tc_analyze
            rows = tc_analyze.analyze(spikes, traces, meta)
            print("\n=== Per-layer sleep-spindle analysis (sigma 10-15 Hz) ===")
            print(tc_analyze._fmt_table(rows))
            md = outdir / f"{tag}_analysis.md"
            tc_analyze._write_md(rows, meta, md)
            print(f"\nWrote {md}")
        except Exception as e:
            print(f"(analysis failed: {e})")

    print(f"\nAll outputs in {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
