#!/usr/bin/env python
"""Histogram a remage LH5 detector output into a PNG plot.

Reads the ragged edep arrays directly.

    python spectrum.py FILE.lh5 [DETECTOR=all] [EMAX_KEV=3000] [BINWIDTH_KEV=5]

DETECTOR "all" (default) aggregates every stp/det* group into one spectrum.

Output (written next to this script, i.e. the project dir):
    {stem}_{det}_spectrum.png   summed-energy histogram, Th-232/U-238 lines tagged in red
"""

import os
import sys

import h5py
import numpy as np

try:
    import hdf5plugin
except ImportError:
    pass

fname = sys.argv[1]
det = sys.argv[2] if len(sys.argv) > 2 else "all"
emax = float(sys.argv[3]) if len(sys.argv) > 3 else 3000
bw = float(sys.argv[4]) if len(sys.argv) > 4 else 5


def det_esum(f, name):
    # edep is a ragged array (one row per event): sum each event's steps
    g = f[f"stp/{name}/edep"]
    fd = g["flattened_data"][:]
    cl = g["cumulative_length"][:]
    starts = np.r_[0, cl[:-1]]
    esum = np.add.reduceat(fd, starts)
    return esum[np.diff(np.r_[0, cl]) > 0]  # drop events with no steps


with h5py.File(fname) as f:
    dets = sorted(k for k in f["stp"] if k.startswith("det")) if det == "all" else [det]
    esum = np.concatenate([det_esum(f, d) for d in dets])

edges = np.arange(0, emax + bw, bw)
counts, _ = np.histogram(esum, bins=edges)
centers = 0.5 * (edges[:-1] + edges[1:])

# write output next to this script (project dir), not next to the input .lh5
outdir = os.path.dirname(os.path.abspath(__file__))
stem = os.path.basename(fname).replace(".lh5", "")
label = det if det != "all" else f"all {len(dets)} detectors"

# --- PNG plot ---
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib not available — skipping PNG")
    sys.exit(0)

# Th-232 + U-238 chain gamma lines to tag in red.
# (energy keV, chain, isotope, label_x keV, label_y as fraction of ymax)
PEAKS = [
    (238, "Th-232", "Pb-212", 175, 1.13),
    (352, "U-238", "Pb-214", 470, 0.86),
    (583, "Th-232", "Tl-208", 640, 1.13),
    (609, "U-238", "Bi-214", 900, 0.98),
    (1120, "U-238", "Bi-214", 1180, 0.55),
    (1764, "U-238", "Bi-214", 1764, 0.32),
    (2614, "Th-232", "Tl-208", 2614, 0.30),
]

pngout = os.path.join(outdir, f"{stem}_{det}_spectrum.png")
fig, ax = plt.subplots(figsize=(11, 6))
ax.stairs(counts, edges, fill=True, color="#2a78d6")
ax.set_xlabel("summed energy per event [keV]")
ax.set_ylabel(f"counts / {bw:g} keV")
ax.set_title(f"{os.path.basename(fname)} — {label} ({len(esum)} entries with hits)")
ax.set_xlim(0, emax)
ymax = counts.max() if counts.max() > 0 else 1
ax.set_ylim(0, ymax * 1.18)

# tag known gamma lines (only those within the plotted range)
for e, chain, iso, lx, lyf in PEAKS:
    if e > emax:
        continue
    win = (centers > e - 15) & (centers < e + 15)  # peak height near the line
    h = counts[win].max() if win.any() else 0
    ax.annotate(
        f"{e} keV\n{iso} ({chain})",
        xy=(e, h),
        xytext=(lx, ymax * lyf),
        color="red",
        fontsize=8,
        fontweight="bold",
        ha="center",
        va="bottom",
        arrowprops=dict(arrowstyle="->", color="red", lw=0.9),
    )

fig.tight_layout()
fig.savefig(pngout, dpi=150)
print("wrote", pngout)
