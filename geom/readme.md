# geom

The geometries, each in its own folder with its viewers and a build check:

| folder | |
| ------ | --- |
| `l1000geom/` | LEGEND-1000: `l1000.gdml` (legend-pygeom-l1000 0.6.0), the geometry every simulation here uses |
| `l200geom/` | LEGEND-200: `l200.gdml` (legend-pygeom-l200 0.11.0), the real array of 101 detectors from the private legend-metadata `62db60c` |

Each holds its GDML (not in git), `*-autopeel.py` (peels it layer by layer as you rotate it), `*-vis.mac` (the
Geant4 viewer), `*-colors.mac`, `toggles/` (shows one material), `check.mac` (the overlap check) and
`autopeel-cache.npz` (the autopeel's meshes).

From `l1000geom/`:

```bash
export LEGEND1000_METADATA=~/Documents/legend-exp/legend1000-metadata
~/venvs/v/bin/legend-pygeom-l1000 -v l1000.gdml
~/venvs/v/bin/python l1000-autopeel.py
~/venvs/v/bin/python l1000-autopeel.py --tube ../../output/tube_geometry.png
remage -i l1000-vis.mac
remage --ignore-warnings -s GDML=l1000.gdml -s SKIP=false -s NPOINTS=1000 -- check.mac 2>&1 | grep "Overlap is detected"
```

From `l200geom/`, the same (`--public-geom` builds a stand-in array from public test data instead):

```bash
export LEGEND_METADATA=~/Documents/legend-metadata
~/venvs/v/bin/legend-pygeom-l200 -v l200.gdml
~/venvs/v/bin/python l200-autopeel.py
remage -i l200-vis.mac
remage --ignore-warnings -s GDML=l200.gdml -s SKIP=false -s NPOINTS=1000 -- check.mac 2>&1 | grep "Overlap is detected"
```

The autopeel's keys: space play/pause · ← → step · ↑ ↓ speed · 0 restart · c start view · r fit · f fly to cursor ·
q quit. l1000's `--holders`, `--tube`, `--uglar` or `--hpge <png>` draws one part alone instead; l200's last layers
sit under the detectors, so look from below. Both autopeels are also in `/Applications/LEGEND.app` (`app/`).
