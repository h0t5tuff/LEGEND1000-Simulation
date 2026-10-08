#  the detector response, with reboost: the steps of a remage run become what LEGEND-1000 would measure, one ROOT file
#  per run for ana/background.C. the recipe and every parameter are the LEGEND-1000 simulation production's
#  (legend-simflow's hit and opt tiers, legend1000-metadata simprod/config, experiment l1000dsg01)
#    ~/venvs/v/bin/python sim/response.py output/tl208.lh5 output/bi214.lh5 ...      -> output/tl208_hit.root ...
#  germanium, per hit: the active energy (FCCD, dead-layer fraction), smeared by the resolution; the A/E classifier from
#  the drift-time map and the current-pulse template. UGLAr, per hit: the photoelectrons the SiPMs see, from the
#  scintillation yield and the optical map

import argparse
import os
import time
import zlib

import awkward as ak
import lh5
import numpy as np
import pint
import pyg4ometry
import pygeomhpges
import pygeomtools
import reboost
import reboost.hpge
import reboost.math
import reboost.spms
import uproot
from dbetto.utils import load_dict
from legendmeta import Legend1000Metadata

# -------------------------------------------------------------------------------
#  1. Configuration:
p = argparse.ArgumentParser(description="detector response of remage runs, with reboost")
p.add_argument("stp", nargs="+", help="remage LH5 outputs (-m: one file per run)")
p.add_argument("--gdml", default="geom/output/l1000.gdml")
p.add_argument("--metadata", default=os.environ.get("LEGEND1000_METADATA", os.path.expanduser("~/Documents/legend1000-metadata")))
p.add_argument("--optmap", default="output/merged_optmap_20260225_063750.lh5", help="UGLAr optical map; its group /all, every SiPM summed, is used")
p.add_argument("--dtmap", default="output/dtmap_V00000A.lh5", help="drift-time maps of the detector template")
p.add_argument("--seed", type=int, default=1)
args = p.parse_args()

EXP = "l1000dsg01"
cfgdir = f"{args.metadata}/simprod/config"
pars = lambda what: load_dict(f"{cfgdir}/pars/{EXP}/geds/{what}/l1000-p01-r%-T%-all-{what}.yaml")  # noqa: E731
hitset = load_dict(f"{cfgdir}/tier/hit/{EXP}/settings.yaml")
optset = load_dict(f"{cfgdir}/tier/opt/{EXP}/settings.yaml")
eres, aoeres, currmod = pars("eresmod"), pars("aoeresmod"), pars("currmod")
dlf = hitset["dead_layer_fraction"]
diodes = Legend1000Metadata(args.metadata, lazy=True).hardware.detectors.germanium.diodes


def per_det(table, det):  # a detector's entry, or the default
    return table.get(det, table["default"])


def fwhm(det, e):  # FWHMLinear: sqrt(a + b E) [keV]
    q = per_det(eres, det)["parameters"]
    return np.sqrt(q["a"] + q["b"] * e)


def aoe_sigma(det, e):  # SigmaFit: sqrt(a + (b / E)^c)
    q = per_det(aoeres, det)["parameters"]
    return np.sqrt(q["a"] + (q["b"] / (e + 1e-99)) ** q["c"])


def smear(x, sigma, rng):  # gaussian, never below zero (legend-simflow's gauss_smear)
    y = np.asarray(x, dtype=float) + rng.standard_normal(len(x)) * np.asarray(sigma, dtype=float)
    return np.where((y <= 0) & (np.asarray(x) >= 0), np.finfo(float).tiny, y)


# -------------------------------------------------------------------------------
#  2. The inputs every run shares: the detectors from the GDML, the drift-time maps, the optical map
t_start = time.time()
reg = pyg4ometry.gdml.Reader(args.gdml).getRegistry()
senstables = pygeomtools.detectors.get_all_senstables(reg)
geds = {n: m for n, m in senstables.items() if m.detector_type == "germanium"}
hpge = {n: pygeomhpges.make_hpge(m.metadata, registry=None, allow_cylindrical_asymmetry=False) for n, m in geds.items()}
dtmaps = reboost.hpge.load_hpge_drift_time_maps(args.dtmap, lh5.ls(args.dtmap)[0], bounds_error=False)
optmap = reboost.spms.load_optmap(args.optmap, "all")
# the map was made in the pygeom-l1000 v0.4.0 argon, 68 mm narrower than this tube's: a step beyond its last radius
# with statistics, at its height, is read there (the detection probability is flat in r up to that edge)
ex, ey, ez = optmap.edges
has = optmap.weights[0] >= 0  # -1: no statistics
cx, cy = (ex[1:] + ex[:-1]) / 2, (ey[1:] + ey[:-1]) / 2
R = np.hypot(*np.meshgrid(cx, cy, indexing="ij"))
r_cov = np.array([R[has[:, :, k]].max() if has[:, :, k].any() else 0.0 for k in range(has.shape[2])]) - 0.005  # m
print(f"inputs: {len(geds)} detectors, drift-time maps {sorted(dtmaps)} deg, optical map {optmap.weights.shape[1:]} bins "
      f"x [{ex[0]}, {ex[-1]}] y [{ey[0]}, {ey[-1]}] z [{ez[0]}, {ez[-1]}] m, read in {time.time() - t_start:.0f} s")
print(f"response: FCCD from the metadata, dead-layer fraction {dlf}, resolution and A/E from simprod {EXP}; "
      f"optical map scaling {optset['optmap_scaling_factor']}, PE resolution {optset['photoelectron_resolution_sigma']}, "
      f"time resolution {optset['time_resolution_in_ns']} ns")
u = pint.UnitRegistry()

# -------------------------------------------------------------------------------
#  3. One run:
for stp in args.stp:
    t0 = time.time()
    # one stream per file, from its name: a file comes out the same whatever else is processed with it. the germanium
    # smearing draws from it; reboost samples the photons and photoelectrons with its own, unseeded generator
    rng = np.random.default_rng([args.seed, zlib.crc32(os.path.basename(stp).encode())])
    out = stp.removesuffix(".lh5") + "_hit.root"
    origins = lh5.read("detector_origins", stp)
    on_disk = set(lh5.ls(stp, "stp/"))
    vtx = lh5.read("vtx", stp).view_as("ak")
    G = {k: [] for k in ("evtid", "det", "edep", "energy", "aoe_class")}
    n_nan = 0
    for det, meta in geds.items():  # the germanium
        if f"stp/{det}" not in on_disk:
            continue
        fccd = diodes[det].characterization.combined_0vbb_analysis.fccd_in_mm.value
        cm = per_det(currmod, det)
        tmpl, times = reboost.hpge.get_current_template(low=-1000, high=4000, step=1, mean_aoe=1, **cm["current_pulse_pars"])
        loc = [origins[det][f].value for f in ("xloc", "yloc", "zloc")] * u.m
        h = lh5.read(f"stp/{det}", stp).view_as("ak", with_units=True)
        dist = reboost.hpge.distance_to_surface(h.xloc, h.yloc, h.zloc, hpge[det], loc, distances_precompute=h.dist_to_surf,
                                                precompute_cutoff=fccd + 1, surface_type="nplus")
        edep_active = h.edep * reboost.math.piecewise_linear_activeness(dist, fccd_in_mm=fccd, dlf=dlf)
        e_true = np.asarray(ak.sum(edep_active, axis=-1), dtype=float)
        energy = smear(e_true, fwhm(det, e_true) / 2.35482, rng)
        dt = reboost.hpge.drift_time_crystal_axes(h.xloc, h.yloc, h.zloc, dtmaps, coord_offset=loc)
        a_max = smear(np.asarray(reboost.hpge.maximum_current(edep_active, dt, template=tmpl, times=times), dtype=float),
                      cm["current_reso"] / cm["mean_aoe"], rng)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):  # hits all in the dead layer: no A/E
            aoe_corr = a_max / energy  # the A/E mean model is a*E + b with a = 0, b = 1: no correction
            aoe_class = np.where(e_true > 0, (aoe_corr - 1) / aoe_sigma(det, e_true), np.nan)
        n_nan += int(np.isnan(aoe_class[e_true > 0]).sum())
        for k, v in (("evtid", h.evtid), ("det", np.full(len(h), meta.uid)), ("edep", ak.sum(h.edep, axis=-1)), ("energy", energy),
                     ("aoe_class", aoe_class)):
            G[k].append(np.asarray(v))
    G = {k: np.concatenate(v) if v else np.zeros(0) for k, v in G.items()}

    L = {k: np.zeros(0) for k in ("evtid", "edep", "pe", "pe_map")}
    note = "no UGLAr steps"
    if "stp/liquid_argon_underground" in on_disk:  # the argon
        a = lh5.read("stp/liquid_argon_underground", stp).view_as("ak")
        x, y, z = a.xloc, a.yloc, a.zloc
        k = np.clip(np.digitize(ak.to_numpy(ak.flatten(z)), ez) - 1, 0, len(r_cov) - 1)
        r = np.hypot(ak.to_numpy(ak.flatten(x)), ak.to_numpy(ak.flatten(y)))
        f = np.where(r > r_cov[k], r_cov[k] / np.maximum(r, 1e-9), 1.0)
        moved = ak.unflatten(f < 1, ak.num(x))
        xc, yc = ak.unflatten(ak.flatten(x) * f, ak.num(x)), ak.unflatten(ak.flatten(y) * f, ak.num(y))
        ph = reboost.spms.emitted_scintillation_photons(a.edep, a.particle, "lar")
        pe = {}
        for key, (px, py) in (("pe", (xc, yc)), ("pe_map", (x, y))):
            n_pe, _, stats = reboost.spms.number_of_detected_photoelectrons(
                px, py, z, ph, optmap, "all", map_scaling=optset["optmap_scaling_factor"],
                max_pes_per_hit=optset["max_pes_per_hit_combined"], return_stats=True)
            tpe = ak.sort(reboost.spms.photoelectron_times(n_pe, a.particle, a.time, "lar"), axis=-1)
            amp = reboost.spms.smear_photoelectrons(tpe, optset["photoelectron_resolution_sigma"], rng=rng)
            _, amp = reboost.spms.cluster_photoelectrons(tpe, amp, optset["time_resolution_in_ns"])
            pe[key] = ak.to_numpy(ak.sum(amp, axis=-1))
            if key == "pe":
                lost = 100 * (stats.energy_oob + stats.energy_no_stats) / max(stats.energy_looped, 1e-12)
        e_all = ak.sum(a.edep)
        note = (f"UGLAr energy read at the map's edge {100 * ak.sum(a.edep[moved]) / e_all:.1f}%, "
                f"outside the map even so {lost:.1f}% (above z {ez[-1]} m: no light)")
        L = {"evtid": ak.to_numpy(a.evtid), "edep": ak.to_numpy(ak.sum(a.edep, axis=-1)), "pe": pe["pe"], "pe_map": pe["pe_map"]}

    trees = {"geds": {"evtid": G["evtid"].astype(np.int64), "det": G["det"].astype(np.int32), "edep": G["edep"].astype(np.float32),
                      "energy": G["energy"].astype(np.float32), "aoe_class": G["aoe_class"].astype(np.float32)},
             "lar": {"evtid": L["evtid"].astype(np.int64), "edep": L["edep"].astype(np.float32), "pe": L["pe"].astype(np.float32),
                     "pe_map": L["pe_map"].astype(np.float32)},
             "vtx": {"evtid": ak.to_numpy(vtx.evtid).astype(np.int64), "xloc": ak.to_numpy(vtx.xloc).astype(np.float32),
                     "yloc": ak.to_numpy(vtx.yloc).astype(np.float32), "zloc": ak.to_numpy(vtx.zloc).astype(np.float32)}}
    with uproot.recreate(out) as f:  # TTrees: uproot's default is now RNTuple
        for name, cols in trees.items():
            f.mktree(name, {k: v.dtype for k, v in cols.items()}).extend(cols)
    print(f"{stp} -> {out}: {len(vtx)} decays, {len(G['evtid'])} germanium hits (A/E undefined for {n_nan}), "
          f"{len(L['evtid'])} UGLAr hits; {note}; {time.time() - t0:.0f} s")
