# geom

The LEGEND-1000 geometry, its viewers and a build check. Run from this folder.

| file | |
| ---- | --- |
| `output/l1000.gdml` | the geometry: legend-pygeom-l1000 0.6.0 at its default `radiogenic` detail, from legend1000-metadata `6cd0209` (not in git) |
| `output/l1000-parts.yaml` | every part of it: material, placements, mass (solid minus daughters, from a 100-slice mesh) |
| `l1000-autopeel.py` | VTK viewer: peels the geometry layer by layer while you rotate it |
| `l1000-vis.mac` | Geant4 viewer: hide and show volumes by hand |
| `l1000-colors.mac` | one colour per volume, for both viewers |
| `toggles/<material>.mac` | shows every volume of one material |
| `check.mac` | does a GDML build, does anything overlap |
| `output/autopeel-cache.npz` | the autopeel's meshed geometry |

```bash
export LEGEND1000_METADATA=~/Documents/legend1000-metadata
~/venvs/v/bin/legend-pygeom-l1000 -v output/l1000.gdml
~/venvs/v/bin/legend-pygeom-l1000 -v --write-manifest output/l1000-parts.yaml
~/venvs/v/bin/python l1000-autopeel.py
remage -i l1000-vis.mac
remage --ignore-warnings -s GDML=output/l1000.gdml -s SKIP=false -s NPOINTS=1000 -- check.mac 2>&1 | grep "Overlap is detected"
```

- **The geometry:** 30 s, the manifest 58 s (7 Oct 2026): 386 logical volumes, 385 parts, 374 233 kg.
  The GDML carries the detector map remage reads with `RegisterDetectorsFromGDML`, so it needs no
  detector macro. Names: the tube `reentrance_tube_copper` (the EFCu, mother of the other two),
  `reentrance_tube_layer_copper_ofhc`, `reentrance_tube_layer_steel_316L`; the argon
  `liquid_argon_atmospheric`, `liquid_argon_underground`; the detectors `V00101Z` ... `V04208Z`
  (string, position); the holders `hpge_string_support_weldment_copper`, placed 1008 times. Both
  coppers are `metal_copper` (8.96 g/cm³), the steel `metal_steel_316L` (8.0 g/cm³).
- **The tube's SS : Cu : EFCu as built:** 619.0 : 669.9 : 231.3 kg, seams at z 920 and 2925 mm.
- **Autopeel:** space play/pause · ← → step · ↑ ↓ speed · 0 restart · c start view · r fit · f fly to cursor ·
  q quit. The first run meshes the GDML (~20 s), later runs start in a second. The peel order is
  `LAYERS`: 18 layers, the 336 detectors stay. Double-clicking `/Applications/l1000.app` runs a standalone
  copy of it, geometry included, that runs on any Mac; `geom/apps/make_app.py l1000` builds it (see `geom/apps/README.md`).
- **Geant4 viewer:** `/vis/geometry/set/visibility <LV> 0 false` hides a volume, `/control/execute
  toggles/pen.mac` shows a material, `/vis/geometry/set/visibility world -1 true` shows everything.
- **Overlap check:** prints nothing when nothing overlaps: none on 7 Oct 2026, in 83 s. `SKIP=true` only
  builds (9 s).
