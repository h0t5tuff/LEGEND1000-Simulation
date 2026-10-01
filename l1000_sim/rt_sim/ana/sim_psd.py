#  step 2: detector response - dead layer and pulse-shape discrimination (A/E), with reboost
#    ~/venvs/v/bin/python ana/sim_psd.py output/tl208.root          # seed 1
#    ~/venvs/v/bin/python ana/sim_psd.py output/tl208.root 7        # another seed
#  writes output/tl208_psd.csv, one row per detector response: evtid,det,energy_keV,aoe_class
#  this is Edgar's build_hit.yaml "geds" group, operation for operation, using reboost's own
#  functions. it reads remage's ROOT output with uproot, so it runs where remage has no HDF5

import csv
import os
import sys

import awkward as ak
import numpy as np
import uproot
from reboost.hpge import psd, utils
from reboost.math import functions, stats
from reboost.units import ureg as u

# ------------------------------------------------------------------------------
#  1. Parameters, from Edgar's reboost/config (read-only mirror):
HERE = os.path.dirname(os.path.abspath(__file__))
REBOOST = os.path.join(HERE, "..", "..", "Edgars_sim", "l1000_simulations", "reboost")
DT_MAP = os.path.join(REBOOST, "drift_time_maps", "drift_time_map.lh5")
REF = "V00000A"                  # every detector uses this one reference set, as in Edgar's config

FCCD_MM, DLF = 1.0, 0.5          # pars.yaml: full-charge-collection depth and dead-layer fraction
PULSE = dict(amax=797, mu=9.78, sigma=60, tail_fraction=0.4,           # pulse_pars.yaml
             tau=300, high_tail_fraction=0.0, high_tau=100.0)
MEAN_AOE, CURRENT_RESO = 0.4, 6
AOE_MU = dict(c=1, m=0)                                                # aoe_class_paras.yaml
AOE_SIGMA = dict(a=2.22383e-06, b=2.44468e01, c=2.27818e00)

# ------------------------------------------------------------------------------
#  2. Detector-independent objects, built once:
load_field = getattr(utils, "load_hpge_rz_field", None) or utils.get_hpge_rz_field  # renamed in reboost 1.x, same arguments
dt_45 = load_field(DT_MAP, REF, "drift_time_045_deg", out_of_bounds_val=0,
                                 method="linear", bounds_error=False)  # slow crystal axis
dt_00 = load_field(DT_MAP, REF, "drift_time_000_deg", out_of_bounds_val=0,
                                 method="linear", bounds_error=False)  # fast crystal axis
template = psd.get_current_template(mean_aoe=1, **PULSE)[0]
times = psd.get_current_template(mean_aoe=MEAN_AOE, **PULSE)[1]


def by_event(evtid, *fields):
    """Sort one detector's steps by event and cut them into one list per event."""
    order = np.argsort(evtid, kind="stable")
    ev = evtid[order]
    starts = np.flatnonzero(np.r_[True, ev[1:] != ev[:-1]])
    counts = np.diff(np.r_[starts, len(ev)])
    return ev[starts], [ak.unflatten(f[order], counts) for f in fields]


def respond(steps, origin, seed):
    """One detector's steps -> active energy and AoE_class per event."""
    ev, (edep, x, y, z, dsurf) = by_event(*steps)
    # dead layer: remage's distance to the nearest surface stands in for the distance to the
    # n+ contact. they differ only near the p+ contact and the passivated groove
    active = functions.piecewise_linear_activeness(dsurf * 1000, fccd_in_mm=FCCD_MM, dlf=DLF)
    e_act = edep * active
    energy = ak.to_numpy(ak.sum(e_act, axis=-1))

    # drift time on both crystal axes, blended by azimuth. positions stay in metres to match
    # the map, so the offset is metres too; remage arrays carry no units, so none are converted
    off = origin * u.m
    t45 = psd.drift_time(x, y, z, dt_45, coord_offset=off)
    t00 = psd.drift_time(x, y, z, dt_00, coord_offset=off)
    phi = np.arctan2(y - origin[1], x - origin[0])
    dtc = t00 + (t45 - t00) * (1 - np.cos(4 * phi)) / 2

    amax = psd.maximum_current(e_act, dtc, template=template, times=times, return_mode="current")
    amax = stats.gaussian_sample(amax, CURRENT_RESO / MEAN_AOE, seed=seed)  # seeded: the same run gives the same PSD
    with np.errstate(divide="ignore", invalid="ignore"):
        aoe = ak.to_numpy(amax) / energy
        aoe_class = (aoe / (AOE_MU["c"] + AOE_MU["m"] * energy) - 1) / np.sqrt(
            AOE_SIGMA["a"] + (AOE_SIGMA["b"] / energy) ** AOE_SIGMA["c"])
    return ev, energy, aoe_class


# ------------------------------------------------------------------------------
#  3. Every germanium detector in the run:
def main(path, seed=1):
    stp = uproot.open(path)["stp"]
    o = stp["detector_origins"].arrays(["name", "xloc_in_m", "yloc_in_m", "zloc_in_m"], library="np")
    origin = {str(n): np.array([x, y, z]) for n, x, y, z in  # rows repeat once per thread
              zip(o["name"], o["xloc_in_m"], o["yloc_in_m"], o["zloc_in_m"])}

    germanium = sorted(k.split(";")[0] for k in stp.keys()
                       if len(k.split(";")[0]) == 5 and k[0] == "V" and k[1:5].isdigit())
    out = os.path.join("output", os.path.basename(path).replace(".root", "") + "_psd.csv")
    os.makedirs("output", exist_ok=True)
    rows = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["evtid", "det", "energy_keV", "aoe_class"])
        for i, det in enumerate(germanium):
            t = stp[det]
            if t.num_entries == 0 or det not in origin:
                continue
            a = t.arrays(["evtid", "edep_in_keV", "xloc_in_m", "yloc_in_m", "zloc_in_m",
                          "dist_to_surf_in_m"], library="np")
            ev, energy, aoe_class = respond(
                (a["evtid"], a["edep_in_keV"], a["xloc_in_m"], a["yloc_in_m"], a["zloc_in_m"],
                 a["dist_to_surf_in_m"]), origin[det], seed * 1000 + i)  # one stream per detector, so their noise is independent
            for i in range(len(ev)):
                if energy[i] > 0:                       # an all-dead-layer response has no A/E
                    w.writerow([int(ev[i]), det, f"{energy[i]:.4f}", f"{aoe_class[i]:.5f}"])
                    rows += 1
    print(f"{path}: {len(germanium)} detectors, {rows} responses, seed {seed} -> {out}")


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit("usage: python ana/sim_psd.py output/<run>.root [seed]")
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) == 3 else 1)
