# LEGEND1000-Simulation

Geant4 simulations, run with [remage](https://github.com/legend-exp/remage), for LEGEND's search for
neutrinoless double-beta decay: what the materials around the germanium detectors (HPGe) contribute to
the background in LEGEND-1000, and how a single detector is characterized at the HADES test stand.

## What is here

| folder | what it does |
| ------ | ------------ |
| [`l1000_sim/rt/`](l1000_sim/rt/) | the re-entrant tube: its Tl-208 and Bi-214 background, and the cheapest steel / copper / electroformed copper (EFCu) tube that stays under budget |
| [`l1000_sim/EFCuHolders/`](l1000_sim/EFCuHolders/) | the EFCu plates holding each detector: the Th-232 and U-238 spectrum across all 336 detectors |
| [`HADES/`](HADES/) | detector characterization with Am-241, Ba-133, Co-60 and Th-228 sources: dead layer, active volume, resolution |
| [`l1000_geometry/`](l1000_geometry/) | an interactive viewer for the LEGEND-1000 geometry: peel it, toggle volumes and materials |
| [`LEGEND200/`](LEGEND200/) | the LEGEND-200 geometry, with the same viewer |
| `l1000_sim/Edgars_sim/`, `l1000_sim/Ralphs_sim/` | colleagues' work the simulations draw on: Edgar's detector response and survival fractions, Ralph's material activities. Kept locally, not in git |

Each simulation folder has its own README: how to run it, what every file does, and its results.

## Geometry

| file | |
| ---- | --- |
| `l1000.gdml`, `l1000-dets.mac` | LEGEND-1000, as built by legend-pygeom-l1000 |
| `KSendcap_l1kGeometry.gdml` | LEGEND-1000 with the KS re-entrant tube endcap, used by `rt` |
| `LEGEND200/l200.gdml`, `l200-dets.mac` | LEGEND-200, as built by legend-pygeom-l200 |

The geometry files are large and not in git: build them with legend-pygeom-l1000 or
legend-pygeom-l200 (which need legend-metadata access) and place them at these paths. `HADES`
builds its own.

## Tools

remage v0.26.0 (Geant4 11.4.2), ROOT 6.40 for the analysis macros, and a Python environment with
the legend-pygeom packages and reboost. Long runs are set up for NERSC (see `rt`).

`rt` and `HADES` keep their source in `sim/` (the simulation) and `ana/` (the analysis of its output),
and write everything they generate to `output/`.

## License

MIT, see [LICENSE](LICENSE).
