# LEGEND200

The LEGEND-200 geometry, its viewers and a build check, the same tools as `geom/`. Run from this folder.

| file | |
| ---- | --- |
| `l200.gdml` | the geometry: legend-pygeom-l200 0.11.0 with its defaults, from legend-metadata `62db60c` (v1.5.5) (not in git) |
| `l200-autopeel.py` | VTK viewer: peels the geometry layer by layer while you rotate it |
| `l200-vis.mac` | Geant4 viewer: hide and show volumes by hand |
| `l200-colors.mac` | one colour per volume, for both viewers |
| `toggles/<material>.mac` | shows every volume of one material (`EnrichedGermanium.mac`: all 101 detectors) |
| `check.mac` | does a GDML build, does anything overlap |
| `output/autopeel-cache.npz` | the autopeel's meshed geometry |

```bash
export LEGEND_METADATA=~/Documents/legend-metadata
~/venvs/v/bin/legend-pygeom-l200 -v l200.gdml
~/venvs/v/bin/python l200-autopeel.py
remage -i l200-vis.mac
remage --ignore-warnings -s GDML=l200.gdml -s SKIP=false -s NPOINTS=1000 -- check.mac 2>&1 | grep "Overlap is detected"
```

- **The geometry:** 15 s (8 Oct 2026). The real array, from the private legend-metadata at the package's
  default timestamp `20230311T235840Z`: 101 HPGe detectors (142.89 kg, 24 enrichments), 58 SiPM
  channels, one volume per fibre; 512 logical volumes, 32 245 placements. The GDML carries the detector
  map: `/RMG/Geometry/RegisterDetectorsFromGDML Germanium` (and `Optical`, `Scintillator` for the SiPMs
  and `liquid_argon`) registers every detector, so no detector macro is needed. Most SiPM channels are
  several volumes (`S054_0`, `S054_1`, ...) sharing the channel's uid. `--public-geom` builds a stand-in
  array from public test data instead, without legend-metadata; a later array needs `metadata_timestamp`
  in a `--config` file.
- **Autopeel:** space play/pause · ← → step · ↑ ↓ speed · 0 restart · c start view · r fit · f fly to cursor ·
  q quit. The first run meshes the GDML (~10 s), later runs start in a second. The peel order is `LAYERS`
  (19 layers, the 101 detectors stay); the last layers (LMFE chips, springs) sit under the detectors, so
  look from below. Double-clicking `/Applications/l200.app` runs it too, its printout in
  `$TMPDIR/l200-autopeel.log`; `geom/apps/make_app.py l200` builds that app (see `geom/apps/README.md`).
- **Geant4 viewer:** `/vis/geometry/set/visibility <LV> 0 false` hides a volume, `/control/execute
  toggles/pen.mac` shows a material, `/vis/geometry/set/visibility world -1 true` shows everything.
- **Overlap check:** prints nothing when nothing overlaps: none on 8 Oct 2026, in 11 s. `SKIP=true` only
  builds (5 s).
