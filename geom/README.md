# l1000_geometry

Two viewers and a build check for the LEGEND-1000 geometry `l1000.gdml` (a link to the repo-level
file). Run everything from this folder.

| file | |
| ---- | --- |
| `l1000-autopeel.py` | VTK window that peels the geometry layer by layer on a timer while you rotate it |
| `l1000-vis.mac` | Geant4 (remage) viewer, for hiding and showing volumes by hand |
| `l1000-colors.mac` | one colour per logical volume (all 2756), read by both viewers |
| `toggles/<material>.mac` | shows every volume of one material (15 materials) |
| `check.mac` | builds the GDML, and optionally scans it for overlaps |
| `output/autopeel-cache.npz` | written by the autopeel: the meshed geometry (22 MB, not in git) |

## Autopeel

```bash
~/venvs/v/bin/python l1000-autopeel.py
```

Every 3 s one layer fades out, outside in, and the camera stays yours. It is a VTK window, not a
Geant4 macro, because remage's Geant4 (11.3, Qt5) does not respond to the mouse while a macro runs.

| key | |
| --- | --- |
| drag, wheel, shift+drag | rotate, zoom, pan |
| space | play / pause (at the end: start again) |
| → / ← | peel the next layer now / put the last one back |
| ↑ / ↓ | peel faster / slower |
| 0 | everything back, start again |
| c / r | start view / fit the camera to what is left (rotation then centres on it) |
| f | fly to the point under the cursor |
| q | quit |

The first run reads the GDML with pyg4ometry, meshes each logical volume once and caches the
surfaces and their placements in `output/` (~25 s). Later runs load the cache and are ready in under
a second; it is rebuilt when `l1000.gdml` changes (delete `output/` to force it). The volumes are
sorted into the layers of `LAYERS` (name regexes; edit it to change the peel order) and drawn as one
actor per layer, colour and hierarchy depth. The timer only changes visibility and opacity. A
daughter sharing a face with its mother is drawn in front of it, so the re-entrant tube's three
sections show without z-fighting. The 1 mm ASICs appear once the PEN plates are gone: pause, look
from below, press `f` on a detector and zoom in.

It prints the layers it found:

```
  #  layer                                                                          LVs  placed
  1  outer cryostat: outer steel tank, skirt and foot                                 3       3
  2  vacuum gap between the two steel tanks                                           1       1
  3  inner cryostat                                                                   1       1
  4  atmospheric LAr (outer bath)                                                     1       1
  5  neutron moderator                                                                1       1
  6  outer reflector: TPB coating                                                     1       1
  7  outer reflector: tetratex lining                                                 1       1
  8  re-entrant tube: EFCu (z < 852 mm), OFHC Cu (852-2857), SS 316L (> 2857)         3       3
  9  underground LAr (inner active bath)                                              1       1
 10  inner reflector: TPB coating                                                     1       1
 11  inner reflector: tetratex lining                                                 1       1
 12  WLS fibres: outer and inner cladding, TPB coating, PS core                       4   48384
 13  SiPM mounts (copper rings, top and bottom)                                       2     252
 14  SiPMs                                                                            1     252
 15  string hardware: rods, weldments, support structure, signal and HV cables     1051    1848
 16  ultem insulators and signal/HV clamps                                         1008    1680
 17  PEN holder plates                                                                2     378
 18  signal ASICs (1 mm silica chips under each detector: look from below)          336     336
  -  the 336 HPGe detectors                                                         336     336

2.78 M triangles in 24 actors, ready in 0.6 s
```

## Geant4 viewer

```bash
remage -i l1000-vis.mac
```

Same colours and start view. In the session:

| command | |
| --- | --- |
| `/vis/geometry/set/visibility <LV> 0 false` | hide one volume (`true` shows it again) |
| `/vis/geometry/set/visibility world -1 false` | hide everything (Geant4 warns "Scene has no extent": harmless) |
| `/control/execute toggles/<material>.mac` | show one material, e.g. `toggles/pen.mac` |
| `/vis/geometry/set/visibility world -1 true` | everything back |

## Overlap check

```bash
remage -q --ignore-warnings -s GDML=l1000.gdml -s SKIP=true -s NPOINTS=1000 -- check.mac
```

```bash
remage --ignore-warnings -s GDML=l1000.gdml -s SKIP=false -s NPOINTS=1000 -- check.mac 2>&1 | grep "Overlap is detected"
```

The first only builds the geometry (7 s). The second also tests 1000 surface points of every volume
(1.5 min) and prints nothing when nothing overlaps, as on 6 Oct 2026. For another GDML, pass
e.g. `-s GDML=../KSendcap_l1kGeometry.gdml`.
