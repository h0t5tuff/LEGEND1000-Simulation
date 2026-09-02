# RT backgrounds — KS-endcap L1000 geometry

| step                             | files                                       |
| -------------------------------- | ------------------------------------------- |
| 1 verify the geometry            | `geo_check.mac`                             |
| 2 verify the reentrance tube     | `rt_verify.C`                               |
| 3 generate events along the tube | `src_events.mac` `src_vertices.C`           |
| 4 run the sims                   | `sim_ge.mac` `sim_optmap.mac` `sim_reach.C` |

`rt_geom.h` is the shared library it reads RT shape from the GDML so no macro hard-codes numbers

```bash
export G=../../KSendcap_l1kGeometry.gdml
```

## 1. Does geometry build, does anything overlap

```bash
remage -q --ignore-warnings -s GDML=$G -s SKIP=true -s NPOINTS=1000 -- geo_check.mac
```

with overlap scan, ~3 min: (for finer scan, increase the points)

```bash
remage --ignore-warnings -s GDML=$G -s SKIP=false -s NPOINTS=1000 -- geo_check.mac 2>&1 | grep "Overlap is detected"
```

## 2. Is RT one clean surface, is the endcap flush

```bash
root -l -b -q rt_verify.C
```

spits `rt_verify.png`

checks on-axis points, duplicate vertices, self-intersections, endcap-vs-wall radius, wall-vs-shell agreement, wall thickness per

## 3. Events along the full RT, in all three sections

`src_events.mac` confines to `reentrancetube` (EFCu), `ofhc_cu` and `ss_316l` together. Fire a small batch and check it before spending real time:

```bash
remage -q --ignore-warnings -t 8 -w -o rt_check.root -s GDML=$G -s NEV=5000 -s SEED=1 -- sim_ge.mac
```

```bash
root -l -b -q 'src_vertices.C("rt_check.root")'
```

Healthy: `EFCu ~12% / OFHC ~44% / SS ~44%`, no empty z-bin, nothing in argon. The
split is uneven by design — sampling is uniform by volume, and the EFCu wall is
1.5 mm against 6 mm of OFHC and SS.

## 4. Run, with or without the optical map

```bash
remage -q --ignore-warnings -t 8 -w -o rt_ge.root -s GDML=$G -s NEV=1000000 -s SEED=1 -- sim_ge.mac
```

```bash
root -l -b -q 'sim_reach.C("rt_ge.root")'
root -l -b -q 'sim_reach.C("rt_ge.root", 24)'   # finer binning
```

`sim_ge.mac` is germanium deposits only — fast, and what you debug with.
`sim_optmap.mac` drops in for the optical-map run: same source, plus the LAr/PEN
deposits and track table, ~9x the size. Neither enables optical physics — the map
replaces photon tracking and is applied afterwards by reboost.

For production, swap the generator lines in `src_events.mac` for the decay chain;
the 2614 keV gamma is a validation placeholder.

## NERSC

Same macros — only `-o` and the launcher change. Locally `.root`, because this
remage reports `HDF5 support: no`. On NERSC `.lh5` for reboost, plus a small
`.root` alongside if you want `src_vertices.C` and `sim_reach.C` on that data.

```bash
remage -q --ignore-warnings -t "$SLURM_CPUS_PER_TASK" -w -o "rt_ge_$SLURM_ARRAY_TASK_ID.lh5" -s GDML=$G -s NEV=1000000 -s SEED="$SLURM_ARRAY_TASK_ID" -- sim_ge.mac
```

Prefix with `apptainer exec --bind <workdir> <image>.sif` if remage is containerised.

## Two traps

`-s` is variadic, so the macro must come after `--`. And remage exits 2 on _any_
warning — every l1000 geometry emits a harmless `SensitiveRegion` one, hence
`--ignore-warnings` on every command above.
