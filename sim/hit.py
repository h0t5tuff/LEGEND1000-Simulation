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
from reboost.optmap import convolve
from dbetto.utils import load_dict
from legendmeta import Legend1000Metadata

# -------------------------------------------------------------------------------
#  1. Configuration: the options, and the response of LEGEND-1000's simulation production (legend1000-metadata simprod)
p = argparse.ArgumentParser(description="detector response of remage runs, with reboost")
p.add_argument("stp", nargs="+", help="remage LH5 outputs (-m: one file per run); the NERSC ones are in output/output_lh5")
p.add_argument("--outdir", default="output/output_root", help="where each <run>_hit.root goes (a link to ~/Downloads/output_root, out of iCloud)")
p.add_argument("--gdml", default="geom/l1000geom/l1000.gdml")
p.add_argument("--metadata", default=os.environ.get("LEGEND1000_METADATA", os.path.expanduser("~/Documents/legend-exp/legend1000-metadata")))
p.add_argument("--optmap", default="output/merged_optmap_20260225_063750.lh5", help="UGLAr optical map; its group /all, every SiPM summed, is used")
p.add_argument("--dtmap", default="output/dtmap_V00000A.lh5", help="drift-time maps of the detector template")
p.add_argument("--seed", type=int, default=1)
p.add_argument("--force", action="store_true", help="rebuild hit files that are already newer than their inputs")
args = p.parse_args()

EXP, CFG = "l1000dsg01", f"{args.metadata}/simprod/config"
pars = lambda what: load_dict(f"{CFG}/pars/{EXP}/geds/{what}/l1000-p01-r%-T%-all-{what}.yaml")  # noqa: E731
eres, aoeres, currmod = pars("eresmod"), pars("aoeresmod"), pars("currmod")
dlf = load_dict(f"{CFG}/tier/hit/{EXP}/settings.yaml")["dead_layer_fraction"]
optset = load_dict(f"{CFG}/tier/opt/{EXP}/settings.yaml")
diodes = Legend1000Metadata(args.metadata, lazy=True).hardware.detectors.germanium.diodes
per_det = lambda table, det: table.get(det, table["default"])  # noqa: E731  a detector's entry, or the default


def smear(x, sigma, rng):  # gaussian, never below zero (legend-simflow's gauss_smear)
    y = np.asarray(x, dtype=float) + rng.standard_normal(len(x)) * np.asarray(sigma, dtype=float)
    return np.where((y <= 0) & (np.asarray(x) >= 0), np.finfo(float).tiny, y)


# -------------------------------------------------------------------------------
#  2. The inputs every run shares: the detectors from the GDML, the drift-time maps, the optical map
t_start = time.time()
senstables = pygeomtools.detectors.get_all_senstables(pyg4ometry.gdml.Reader(args.gdml).getRegistry())
geds = {n: m for n, m in senstables.items() if m.detector_type == "germanium"}
hpge = {n: pygeomhpges.make_hpge(m.metadata, registry=None, allow_cylindrical_asymmetry=False) for n, m in geds.items()}
dtmaps = reboost.hpge.load_hpge_drift_time_maps(args.dtmap, lh5.ls(args.dtmap)[0], bounds_error=False)
optmap = reboost.spms.load_optmap(args.optmap, "all")
SCINT = convolve._get_scint_params("lar")  # what reboost.spms's processors use; called below with our own generator
u = pint.UnitRegistry()
# the map was made in the pygeom-l1000 v0.4.0 argon, 68 mm narrower than this tube's: a step beyond its last radius with
# statistics, at its height, is read there (the detection probability is flat in r up to that edge)
ex, ey, ez = optmap.edges
has = optmap.weights[0] >= 0  # -1: no statistics
R = np.hypot(*np.meshgrid((ex[1:] + ex[:-1]) / 2, (ey[1:] + ey[:-1]) / 2, indexing="ij"))
r_cov = np.array([R[has[:, :, k]].max() if has[:, :, k].any() else 0.0 for k in range(has.shape[2])]) - 0.005  # m
print(f"inputs: {len(geds)} detectors, drift-time maps {sorted(dtmaps)} deg, optical map {optmap.weights.shape[1:]} bins "
      f"x [{ex[0]}, {ex[-1]}] y [{ey[0]}, {ey[-1]}] z [{ez[0]}, {ez[-1]}] m, read in {time.time() - t_start:.0f} s")
print(f"response: FCCD from the metadata, dead-layer fraction {dlf}, resolution and A/E from simprod {EXP}; "
      f"optical map scaling {optset['optmap_scaling_factor']}, PE resolution {optset['photoelectron_resolution_sigma']}, "
      f"time resolution {optset['time_resolution_in_ns']} ns")


# -------------------------------------------------------------------------------
#  3. The germanium: per detector, the energy through the dead layer, smeared by the resolution, and the A/E classifier
def germanium(stp, on_disk, rng):
    origins = lh5.read("detector_origins", stp)
    cols, n_nan = [], 0
    for det, meta in geds.items():
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
        q = per_det(eres, det)["parameters"]  # FWHMLinear: sqrt(a + b E) [keV]
        energy = smear(e_true, np.sqrt(q["a"] + q["b"] * e_true) / 2.35482, rng)
        dt = reboost.hpge.drift_time_crystal_axes(h.xloc, h.yloc, h.zloc, dtmaps, coord_offset=loc)
        a_max = smear(np.asarray(reboost.hpge.maximum_current(edep_active, dt, template=tmpl, times=times), dtype=float),
                      cm["current_reso"] / cm["mean_aoe"], rng)
        q = per_det(aoeres, det)["parameters"]  # SigmaFit: sqrt(a + (b / E)^c)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):  # hits all in the dead layer: no A/E
            # the A/E mean model is a*E + b with a = 0, b = 1: A/E needs no correction
            aoe_class = np.where(e_true > 0, (a_max / energy - 1) / np.sqrt(q["a"] + (q["b"] / (e_true + 1e-99)) ** q["c"]), np.nan)
        n_nan += int(np.isnan(aoe_class[e_true > 0]).sum())
        cols.append({"evtid": h.evtid, "det": np.full(len(h), meta.uid), "edep": ak.sum(h.edep, axis=-1), "energy": energy,
                     "aoe_class": aoe_class})
    G = {k: np.concatenate([np.asarray(c[k]) for c in cols]) if cols else np.zeros(0) for k in ("evtid", "det", "edep", "energy", "aoe_class")}
    return G, f"A/E undefined for {n_nan}"


# -------------------------------------------------------------------------------
#  4. The argon: each step's scintillation photons, then the photoelectrons the SiPMs see through the optical map, twice:
#     pe, with the argon beyond the map's edge read at the edge; pe_map, off the map alone
def argon(stp, on_disk, rng):
    if "stp/liquid_argon_underground" not in on_disk:
        return {k: np.zeros(0) for k in ("evtid", "edep", "pe", "pe_map")}, "no UGLAr steps"
    a = lh5.read("stp/liquid_argon_underground", stp).view_as("ak")
    x, y, z = a.xloc, a.yloc, a.zloc
    k = np.clip(np.digitize(ak.to_numpy(ak.flatten(z)), ez) - 1, 0, len(r_cov) - 1)
    r = np.hypot(ak.to_numpy(ak.flatten(x)), ak.to_numpy(ak.flatten(y)))
    f = np.where(r > r_cov[k], r_cov[k] / np.maximum(r, 1e-9), 1.0)  # pulls a step beyond the edge back onto it
    xc, yc = ak.unflatten(ak.flatten(x) * f, ak.num(x)), ak.unflatten(ak.flatten(y) * f, ak.num(y))
    ph = convolve.iterate_stepwise_depositions_scintillate(ak.Array({"edep": a.edep, "particle": a.particle}), SCINT, rng=rng)
    L = {"evtid": ak.to_numpy(a.evtid), "edep": ak.to_numpy(ak.sum(a.edep, axis=-1))}
    for key, (px, py) in (("pe", (xc, yc)), ("pe_map", (x, y))):
        n_pe, _, stats = convolve.iterate_stepwise_depositions_numdet(
            ak.Array({"xloc": px, "yloc": py, "zloc": z, "num_scint_ph": ph}), optmap, "all", map_scaling=optset["optmap_scaling_factor"],
            max_pes_per_hit=optset["max_pes_per_hit_combined"], rng=rng, return_stats=True)
        tpe = ak.sort(convolve.iterate_stepwise_depositions_times(ak.Array({"num_det_ph": n_pe, "particle": a.particle, "time": a.time}),
                                                                    SCINT, rng=rng), axis=-1)
        amp = reboost.spms.smear_photoelectrons(tpe, optset["photoelectron_resolution_sigma"], rng=rng)
        L[key] = ak.to_numpy(ak.sum(reboost.spms.cluster_photoelectrons(tpe, amp, optset["time_resolution_in_ns"])[1], axis=-1))
        if key == "pe":
            lost = 100 * (stats.energy_oob + stats.energy_no_stats) / max(stats.energy_looped, 1e-12)
    moved = ak.unflatten(f < 1, ak.num(x))
    return L, (f"UGLAr energy read at the map's edge {100 * ak.sum(a.edep[moved]) / ak.sum(a.edep):.1f}%, "
               f"outside the map even so {lost:.1f}% (above z {ez[-1]} m: no light)")


# -------------------------------------------------------------------------------
#  5. One run -> output/output_root/<run>_hit.root: TTrees geds, lar, vtx. One random stream per file, from its name, so a
#     file comes out the same whatever else is processed with it; the argon calls reboost.optmap.convolve directly because
#     reboost.spms's processors draw from an unseeded generator
for stp in args.stp:
    t0 = time.time()
    out = os.path.join(args.outdir, os.path.basename(stp).removesuffix(".lh5") + "_hit.root")
    if not args.force and os.path.exists(out) and os.path.getmtime(out) > max(os.path.getmtime(f) for f in (stp, args.optmap, args.dtmap)):
        print(f"{stp} -> {out}: up to date, kept")  # the same seeds give the same hits: nothing to rebuild
        continue
    rng = np.random.default_rng([args.seed, zlib.crc32(os.path.basename(stp).encode())])
    on_disk = set(lh5.ls(stp, "stp/"))
    G, g_note = germanium(stp, on_disk, rng)
    L, l_note = argon(stp, on_disk, rng)
    vtx = lh5.read("vtx", stp).view_as("ak")
    V = {k: ak.to_numpy(vtx[k]) for k in ("evtid", "xloc", "yloc", "zloc")}
    dtype = lambda k: np.int64 if k == "evtid" else (np.int32 if k == "det" else np.float32)  # noqa: E731
    with uproot.recreate(out) as f:  # TTrees: uproot's default is now RNTuple
        for name, cols in (("geds", G), ("lar", L), ("vtx", V)):
            cols = {k: np.asarray(v).astype(dtype(k)) for k, v in cols.items()}
            f.mktree(name, {k: v.dtype for k, v in cols.items()}).extend(cols)
    print(f"{stp} -> {out}: {len(V['evtid'])} decays, {len(G['evtid'])} germanium hits ({g_note}), {len(L['evtid'])} UGLAr hits; "
          f"{l_note}; {time.time() - t0:.0f} s")
