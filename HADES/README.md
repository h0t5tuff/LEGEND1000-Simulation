# HPGe characterization at HADES

One HPGe in a vacuum cryostat with a calibration source on top: the stand LEGEND characterizes its
detectors in, built by [legend-pygeom-hades](https://legend-pygeom-hades.readthedocs.io) and simulated
with remage. **The unknown is the detector, not the background**: each source is a known input, and
what the crystal makes of it measures its dead layer (FCCD, the full-charge-collection depth), its
active volume and its resolution.

| source | stand | measures | observable | simulated as |
| ------ | ----- | -------- | ---------- | ------------ |
| Am-241 | `am_HS1` | dead layer | 59.5 keV peak per decay | the 59.5 keV line, aimed (the source is collimated) |
| Ba-133 | `ba_HS4` | dead layer | R = (79.6+81.0 keV) / 356.0 keV: no activity needed | decays |
| Co-60  | `co_HS5` | active volume | 1332.5 keV peak per decay | decays |
| Th-228 | `th_HS2` | resolution, pulse-shape peaks | 2614.5 keV peak per decay; lines from 239 to 2615 keV, DEP, SEP | decays, the whole chain |

**Method.** `sim/` fires the source and records the energy every step deposits in the crystal, with
its distance to the surface. `ana/` turns that into what the detector reads out: it removes the dead
layer, smears with the resolution and counts the peaks as on a measured spectrum, over a scan of
FCCDs. The dead layer is applied in `ana/`, not `sim/`, so one remage run serves every FCCD, and a
measured observable reads back as the detector's FCCD.

## Files

| source | |
| ------ | --- |
| [`sim/<src>.yaml`](#simsrcyaml) | the stand, per source `<src>` = `am241`, `ba133`, `co60`, `th228`: detector, source, source position |
| [`sim/run.mac`](#simrunmac) | remage: decay the nuclide (or fire one gamma line into a cone) in the source, record every step in the crystal |
| [`ana/spectrum.C`](#anaspectrumc) | steps → detector events → dead layer → resolution → peaks → FCCD scan → a measurement read back as an FCCD |

| output, per source | made by | holds | size |
| ------------------ | ------- | ----- | ---- |
| [`output/<src>.gdml`](#outputsrcgdml) | legend-pygeom-hades | the stand | 16-26 kB |
| [`output/<src>.root`](#outputsrcroot) | remage | every step in the crystal, every primary | 155-281 MB |
| [`output/<src>_spectrum.png`](#outputsrc_spectrumpng) | `spectrum.C` | the spectrum and the FCCD scan | ~37 kB |

`.gitignore` keeps `output/` out of git, except the PNGs. No Python of our own.

## Run

From `HADES/`, with legend-pygeom-hades 0.2.1 (`~/venvs/v/bin/pip install legend-pygeom-hades`),
remage v0.26.0 (Geant4 11.4.2) and ROOT 6.40:

```bash
mkdir -p output
~/venvs/v/bin/legend-pygeom-hades --public-geom --config sim/am241.yaml output/am241.gdml
~/venvs/v/bin/legend-pygeom-hades --public-geom --config sim/ba133.yaml output/ba133.gdml
~/venvs/v/bin/legend-pygeom-hades --public-geom --config sim/co60.yaml output/co60.gdml
~/venvs/v/bin/legend-pygeom-hades --public-geom --config sim/th228.yaml output/th228.gdml
remage -q --ignore-warnings -t 8 -w -o output/am241.root -s GDML=output/am241.gdml -s PARTICLE=gamma -s Z=95 -s A=241 -s E=59.5409 -s THETA=10 -s NEV=4000000 -s SEED=1 -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/ba133.root -s GDML=output/ba133.gdml -s PARTICLE=ion -s Z=56 -s A=133 -s E=0 -s THETA=180 -s NEV=2000000 -s SEED=4 -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/co60.root -s GDML=output/co60.gdml -s PARTICLE=ion -s Z=27 -s A=60 -s E=0 -s THETA=180 -s NEV=2000000 -s SEED=3 -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/th228.root -s GDML=output/th228.gdml -s PARTICLE=ion -s Z=90 -s A=228 -s E=0 -s THETA=180 -s NEV=5000000 -s SEED=2 -- sim/run.mac
root -l -b -q 'ana/spectrum.C("output/am241.root")'
root -l -b -q 'ana/spectrum.C("output/ba133.root")'
root -l -b -q 'ana/spectrum.C("output/co60.root")'
root -l -b -q 'ana/spectrum.C("output/th228.root")'
```

| step | Am-241 | Ba-133 | Co-60 | Th-228 |
| ---- | ------ | ------ | ----- | ------ |
| legend-pygeom-hades | ~12 s | ~12 s | ~12 s | ~12 s |
| remage, `-t 8` | 4 × 10⁶ gammas, 31 s | 2 × 10⁶ decays, 86 s | 2 × 10⁶ decays, 114 s | 5 × 10⁶ decays, 261 s |
| `spectrum.C` | 3 s | 4 s | 4 s | 6 s |

- **More statistics:** raise `NEV`, or add a run with its own `SEED` and file name.
- **What they print:** legend-pygeom-hades one warning, `CONSTRUCTING GEOMETRY FROM PUBLIC DATA ONLY`
  (`--public-geom`). remage one warning, `Could not find G4Region 'SensitiveRegion'` (Geant4's
  default production cuts apply everywhere), then its `G4GDML: Reading ...` lines and `Stripping off
  GDML names of materials, solids and volumes ...`. `spectrum.C` prints its
  [report](#what-spectrumc-prints).

## Source code

### `sim/<src>.yaml`

The stand, read by legend-pygeom-hades. The four files differ only in their first comment line and
`measurement`:

| key | value | meaning |
| --- | ----- | ------- |
| `detector` | `V99000A` | the diode. With `--public-geom`, a dummy from the package: `V` an ICPC (`B99000A` a BEGe) |
| `measurement` | `am_HS1_top_dlt`, `ba_HS4_top_dlt`, `co_HS5_top_dlt`, `th_HS2_top_psa` | `<source>_<HSX>_<position>_<ID>`: builds the source, its holder, and the lead castle (not for `am_HS1`). Sources `am_HS1`, `th_HS2`, `co_HS5`, `ba_HS4`, `am_HS6`; position `top`, or `lat` for `am_HS1` and `th_HS2`; the ID labels the campaign |
| `source_position` | φ 0°, r 0 mm, z 38 mm | where the source sits on the cryostat, in cylindrical coordinates |
| `run` | `run0001` | the campaign run: it changes the castle table only for V02160A's early `th_HS2_lat_psa` runs |
| `daq_settings.flashcam.card_interface` | `efb2` | required by the package; the castle table comes from the detector name (table 1 here) |

### `sim/run.mac`

`remage [flags] -s KEY=VALUE ... -- sim/run.mac`: each `-s` fills a `{KEY}` alias. Flags: `-q` only
warnings and errors, `--ignore-warnings` a warning does not fail the exit code, `-t 8` threads, `-w`
overwrite, `-o` the output file.

| alias | decay mode | line mode | meaning |
| ----- | ---------- | --------- | ------- |
| `GDML` | `output/<src>.gdml` | `output/am241.gdml` | the stand |
| `PARTICLE` | `ion` | `gamma` | what is fired |
| `Z`, `A` | the nuclide | the nuclide | `/gps/ion` needs them; a gamma then replaces the nucleus |
| `E` | `0` | the line [keV] | kinetic energy: a nucleus starts at rest |
| `THETA` | `180` | the cone [deg] | emission half-angle from −z, where the crystal is; a decay ignores it |
| `NEV` | decays | gammas | primaries |
| `SEED` | | | the random seed, one per run |

1. **Before init:** load the GDML; register the diode from the map legend-pygeom-hades writes into it
   (uid 1) and name its step table after it (`stp/V99000A`); skip the overlap check (pyg4ometry did it).
2. **Output:** germanium steps that deposit energy, in single precision, with `dist_to_surf`; and the
   primaries (`StorePrimaryParticleInformation`), from which `spectrum.C` reads the source.
3. **Decay physics:** `DaughterNucleusMaxLifetime 1 y`: a daughter living longer ends the decay, so
   Am-241 stops at Np-237 while the Th-228 chain (members up to 3.6 d) runs to Pb-208.
4. **Source:** GPS confined to `Source_PV`, the active disc (not `source_pv`, the holder around it);
   isotropic up to `THETA`; `beamOn {NEV}`.

- **Decays decay:** Geant4 11.4 ignores primaries living longer than a year; remage v0.26 raises that
  to 3 × 10¹⁰ years, so Am-241 (432 y), Ba-133 (10.5 y), Co-60 (5.3 y) and Th-228 (1.9 y) decay.
- **The chain:** remage sets the primary's decay time to 0, but each later member keeps its real
  delay: one Th-228 event spans up to 70 days (6 × 10¹⁵ ns). `spectrum.C` splits it.
- **Line mode for Am-241:** the `am_HS1` channel (Ø1.0 × 25.6 mm, a 2.24° half-angle) passes
  3.8 × 10⁻⁴ of the gammas, so decays would waste 99.96% of the time. The cone is 10°, not 3°:
  [the check](#line-mode-against-decays).

### `ana/spectrum.C`

```
spectrum.C(file, fccd = 1.0, measured = -1, error = 0, dlf = 0.5, seed = 1)
  file       output/<src>.root. the source is read from it, never named
  fccd       the FCCD [mm] the response and the spectrum are given at
  measured   a measured observable and its error, read back as an FCCD. used only if > 0
  dlf        the fraction of the FCCD that is fully dead
  seed       the resolution smearing
```

| section | |
| ------- | --- |
| 1. configuration | `cfg::` below. `activeness(d, fccd, dlf)`: 0 for d ≤ dlf × fccd, linear up to 1 at d = fccd. `scanTree()` reads named columns as doubles, whatever their stored type |
| 2. the source | from `stp/particles`. A nucleus: Z and A pick the source in `cfg::sources`, counts are per decay. A gamma: its energy picks the line (within 0.05 keV), the widest direction from −z gives the cone θ, and per decay = per gamma × (1 − cos θ)/2 × the line's emission probability; a gamma matching no line is counted per gamma. ERROR and stop: no `stp/particles`, a nuclide not in `cfg::sources`, any other primary, no germanium table |
| 3. events | the germanium table, sorted by event and time. A detector event is the steps within 10 µs of its first, which splits the Th-228 chain. Per event: the deposited energy, and the active energy Σ edep × activeness at every FCCD of the scan |
| 4. peaks | one Gaussian draw per event, shared by every FCCD so the scan moves smoothly. Per peak: counts in [e1 − 3σ, e2 + 3σ] minus two side bands of half its width, error √(window + side bands) |
| 5. response | at `fccd`: the energy lost to the dead layer; per peak net counts, per decay, and the FWHM from a χ² fit (a Gaussian on an erfc step, ±4σ, peaks over 300 net counts). When the fitted peaks span a factor 3 in energy: FWHM² = A + B E |
| 6. FCCD scan | the observable at every FCCD. With `measured`: linear interpolation in its interval; error = (measurement ⊕ MC statistics) / local slope |
| 7. draw | `output/<src>_spectrum.png` |

| `cfg::` | value | |
| ------- | ----- | --- |
| `resoA`, `resoB` | 0.865 keV², 0.00225 keV | FWHM(E) = √(A + B E): 1.0 keV at 60 keV, 2.6 keV at 2615 keV. **A placeholder** for the detector's calibration |
| `window_ns` | 10⁴ ns | the detector event window |
| `scan` | 0, 0.25, ..., 2.0, 2.5, 3.0 mm | the FCCDs scanned; `fccd` is added when off the grid |
| `minFit` | 300 | net counts a peak needs before its width is fitted |
| `sources` | below | each source's peaks and observable |

`cfg::sources`: a peak is a window from `e1` to `e2` keV (a doublet when they differ), with its
gammas per decay (used only by line mode; 0 for doublets, sums and escape peaks); `f` = width fitted.
Bold: the observable.

| source | Z, A | peaks [keV] (gammas per decay) | observable |
| ------ | ---- | ------------------------------ | ---------- |
| Am-241 | 95, 241 | **59.5409** (0.3592, f); 26.3446 (0.0231) | 59.5 per decay |
| Ba-133 | 56, 133 | **79.6142-80.9979**; 276.3989 (0.0716, f); 302.8508 (0.1834, f); **356.0129** (0.6205, f); 383.8485 (0.0894, f) | R = (79.6+81.0) / 356.0 |
| Co-60 | 27, 60 | 1173.228 (0.9985, f); **1332.492** (0.9998, f); 2505.720 sum | 1332.5 per decay |
| Th-228 | 90, 228 | 238.632-240.986 (Pb-212, Ra-224); 583.187 (0.3045, f); 727.330 (0.0667, f); 860.557 (0.0448, f); 1592.511 DEP; 1620.50 (0.0147, f); 2103.511 SEP; **2614.511** (0.3585, f) | 2614.5 per decay |

## Outputs

### `output/<src>.gdml`

In all four: a vacuum world; a lab of air, turned 180° so z points up; the aluminium cryostat
`cryo_pv_1` and its vacuum `cavity_pv`, holding the crystal `V99000A` (germanium enriched to 75%
Ge-76), its aluminium holder `holder_pv_1` and polyethylene wrap `wrap_pv_1`. The source is
`Source_PV` in every one, which is why one macro serves all four:

| source | the source, as built | castle |
| ------ | -------------------- | ------ |
| Am-241 | `Source_PV`, a SiO₂ disc Ø1 × 2 mm, in `Source_Encapsulated_PV`, a steel capsule Ø2 × 10 mm, in `source_pv_1`, a 30 × 30 × 65 mm copper collimator with an air channel `Collimator_Beam_PV` Ø1.0 × 25.6 mm | none |
| Ba-133 | `Source_PV`, a polyethylene disc Ø5 × 0.1 mm, in `Source_foil_PV`, a 0.5 mm polyethylene foil Ø26 mm, in `Source_Alring_PV`, an aluminium ring Ø26-30 × 3 mm; acrylic holder `source_holder_pv_1` | aluminium plate `plate_pv_1`, lead castle `castle_pv_1` |
| Co-60 | as Ba-133, with the foil and the ring's opening Ø20 mm | the same |
| Th-228 | `Source_PV`, a SiO₂ disc Ø1 × 1 mm, and `Source_Epoxy_PV`, an epoxy plug Ø1.6 × 2.2 mm, in `Source_Encapsulated_PV`, a steel capsule Ø2 × 7 mm, in `source_pv_1`, a copper block; lead plates `th_plate_pv_1`; acrylic holder `source_holder_pv_1` | the same |

### `output/<src>.root`

Under `stp/`:

| tree | one row per | columns | read by `spectrum.C` |
| ---- | ----------- | ------- | -------------------- |
| `V99000A` | step in the crystal | `evtid`, `particle` (PDG), `edep_in_keV`, `time_in_ns`, `xloc_in_m`, `yloc_in_m`, `zloc_in_m`, `dist_to_surf_in_m` (to the nearest crystal surface) | `evtid`, `edep_in_keV`, `dist_to_surf_in_m`, `time_in_ns` |
| `particles` | primary | `evtid`, `vertexid`, `particle` (PDG), `px_in_MeV`, `py_in_MeV`, `pz_in_MeV`, `ekin_in_MeV` | `particle`, `px_in_MeV`, `py_in_MeV`, `pz_in_MeV`, `ekin_in_MeV` |
| `vtx` | event | `evtid`, `time_in_ns`, `xloc_in_m`, `yloc_in_m`, `zloc_in_m`, `n_part` | no |
| `detector_origins` | thread (8) | `name`, `xloc_in_m`, `yloc_in_m`, `zloc_in_m`: `V99000A` at z = −0.068 m | no |

Floats are single precision; `time_in_ns` is double, which the Th-228 chain needs.

| rows | Am-241 | Ba-133 | Co-60 | Th-228 |
| ---- | ------ | ------ | ----- | ------ |
| `V99000A` | 434 474 | 3 498 645 | 2 881 374 | 4 288 590 |
| `particles`, `vtx` | 4 000 000 | 2 000 000 | 2 000 000 | 5 000 000 |

### What `spectrum.C` prints

The source and its statistics; at the requested FCCD, the energy lost to the dead layer and each
peak's window, net counts, counts per decay and fitted FWHM (the injected one in parentheses); the
resolution fit, when the peaks allow one; the FCCD scan, `<- drawn` at the requested FCCD.

Am-241:

```
output/am241.root: Am-241 59.54 keV line in a 10.00 deg cone: 4000000 gammas = 5.27e+08 in 4 pi = 1.47e+09 decays
  434474 steps -> 33663 detector events (10 us window)

FCCD 1.00 mm, DLF 0.50, seed 1: energy lost to the dead layer 54.35 %
  peak [keV]   window [keV]              net counts              per decay   FWHM [keV]: fitted (injected)
  59.5           58.27 - 60.81        11271 +- 109     7.688e-06 +- 7.4e-08    1.018 +- 0.009 (0.999)

FCCD scan, DLF 0.50: 59.5 per decay against the dead layer
  FCCD [mm]            59.5 per decay
  0.00          2.1992e-05 +- 1.2e-07 
  0.25            1.66e-05 +- 1.1e-07 
  0.50          1.2787e-05 +- 9.5e-08 
  0.75          9.9189e-06 +- 8.4e-08 
  1.00          7.6883e-06 +- 7.4e-08    <- drawn
  1.25          5.9755e-06 +- 6.6e-08 
  1.50          4.6283e-06 +- 5.8e-08 
  1.75          3.5628e-06 +- 5.1e-08 
  2.00          2.7019e-06 +- 4.5e-08 
  2.50          1.6132e-06 +- 3.5e-08 
  3.00          9.7067e-07 +- 2.7e-08 
wrote output/am241_spectrum.png
```

Ba-133:

```
output/ba133.root: Ba-133: 2000000 decays
  3498645 steps -> 320793 detector events (10 us window)

FCCD 1.00 mm, DLF 0.50, seed 1: energy lost to the dead layer 20.06 %
  peak [keV]   window [keV]              net counts              per decay   FWHM [keV]: fitted (injected)
  79.6+81.0      78.31 - 82.30        28646 +- 179     1.432e-02 +- 9.0e-05    -
  276.4         274.85 - 277.95        4767 +- 74      2.384e-03 +- 3.7e-05    1.198 +- 0.017 (1.219)
  302.9         301.27 - 304.43       11918 +- 112     5.959e-03 +- 5.6e-05    1.239 +- 0.010 (1.244)
  356.0         354.37 - 357.66       35857 +- 191     1.793e-02 +- 9.5e-05    1.288 +- 0.006 (1.291)
  383.8         382.17 - 385.52        5213 +- 73      2.606e-03 +- 3.6e-05    1.331 +- 0.015 (1.315)

FCCD scan, DLF 0.50: R = (79.6+81.0) / (356.0) against the dead layer
  FCCD [mm]      79.6+81.0 /decay         356.0 /decay                        R
  0.00                 2.4653e-02           2.0324e-02         1.213 +- 0.0083  
  0.25                 2.1729e-02           2.0088e-02        1.0817 +- 0.0077  
  0.50                 1.8981e-02           1.9471e-02       0.97483 +- 0.0072  
  0.75                 1.6494e-02           1.8708e-02       0.88168 +- 0.0069  
  1.00                 1.4323e-02           1.7928e-02        0.7989 +- 0.0066     <- drawn
  1.25                 1.2481e-02           1.7167e-02       0.72701 +- 0.0063  
  1.50                 1.0825e-02           1.6441e-02        0.6584 +- 0.0061  
  1.75                 9.4055e-03           1.5709e-02       0.59875 +- 0.0058  
  2.00                 8.1905e-03           1.5062e-02       0.54379 +- 0.0056  
  2.50                 6.1895e-03           1.3714e-02       0.45133 +- 0.0053  
  3.00                 4.6780e-03           1.2460e-02       0.37544 +- 0.005   
wrote output/ba133_spectrum.png
```

Co-60:

```
output/co60.root: Co-60: 2000000 decays
  2881374 steps -> 215576 detector events (10 us window)

FCCD 1.00 mm, DLF 0.50, seed 1: energy lost to the dead layer 7.73 %
  peak [keV]   window [keV]              net counts              per decay   FWHM [keV]: fitted (injected)
  1173.2       1170.84 - 1175.61      22109 +- 151     1.105e-02 +- 7.6e-05    1.864 +- 0.011 (1.872)
  1332.5       1329.99 - 1335.00      19633 +- 141     9.817e-03 +- 7.1e-05    1.977 +- 0.012 (1.965)
  2505.7 sum   2502.47 - 2508.97        258 +- 16      1.290e-04 +- 8.1e-06    -

FCCD scan, DLF 0.50: 1332.5 per decay against the dead layer
  FCCD [mm]          1332.5 per decay
  0.00             0.01184 +- 7.7e-05 
  0.25            0.011293 +- 7.6e-05 
  0.50            0.010774 +- 7.4e-05 
  0.75              0.0103 +- 7.2e-05 
  1.00           0.0098165 +- 7.1e-05    <- drawn
  1.25           0.0093405 +- 6.9e-05 
  1.50           0.0088875 +- 6.7e-05 
  1.75             0.00849 +- 6.6e-05 
  2.00           0.0080805 +- 6.4e-05 
  2.50            0.007314 +- 6.1e-05 
  3.00           0.0065895 +- 5.8e-05 
wrote output/co60_spectrum.png
```

Th-228:

```
output/th228.root: Th-228: 5000000 decays
  4288590 steps -> 397453 detector events (10 us window)

FCCD 1.00 mm, DLF 0.50, seed 1: energy lost to the dead layer 8.49 %
  peak [keV]   window [keV]              net counts              per decay   FWHM [keV]: fitted (injected)
  238.6+241.0   237.12 - 242.50       26586 +- 179     5.317e-03 +- 3.6e-05    -
  583.2         581.31 - 585.07       17770 +- 136     3.554e-03 +- 2.7e-05    1.485 +- 0.010 (1.476)
  727.3         725.32 - 729.34        3633 +- 63      7.266e-04 +- 1.3e-05    1.528 +- 0.026 (1.582)
  860.6         858.42 - 862.69        2185 +- 49      4.370e-04 +- 9.8e-06    1.668 +- 0.035 (1.674)
  1592.5 DEP   1589.82 - 1595.20        477 +- 26      9.540e-05 +- 5.3e-06    -
  1620.5       1617.79 - 1623.21        493 +- 26      9.860e-05 +- 5.2e-06    2.160 +- 0.126 (2.124)
  2103.5 SEP   2100.50 - 2106.52       1043 +- 37      2.086e-04 +- 7.4e-06    -
  2614.5       2611.20 - 2617.82       7390 +- 87      1.478e-03 +- 1.7e-05    2.607 +- 0.025 (2.598)
  resolution from 5 peaks, FWHM^2 = A + B E: A 0.867 +- 0.052 keV^2 (injected 0.865), B 2.255e-03 +- 6.5e-05 keV (injected 2.250e-03)

FCCD scan, DLF 0.50: 2614.5 per decay against the dead layer
  FCCD [mm]          2614.5 per decay
  0.00           0.0017866 +- 1.9e-05 
  0.25           0.0016972 +- 1.9e-05 
  0.50           0.0016208 +- 1.8e-05 
  0.75           0.0015504 +- 1.8e-05 
  1.00            0.001478 +- 1.7e-05    <- drawn
  1.25            0.001402 +- 1.7e-05 
  1.50            0.001335 +- 1.6e-05 
  1.75            0.001271 +- 1.6e-05 
  2.00           0.0012036 +- 1.6e-05 
  2.50           0.0010858 +- 1.5e-05 
  3.00            0.000977 +- 1.4e-05 
wrote output/th228_spectrum.png
```

With a measurement, the report ends with one more line. Here R at FCCD 0.80 mm (0.86712, from
`ana/spectrum.C("output/ba133.root", 0.8)`) is fed back ± 0.01 with
`ana/spectrum.C("output/ba133.root", 1.0, 0.86712, 0.01)`:

```
measured 0.86712 +- 0.01  ->  FCCD 0.794 +- 0.037 mm   (measurement 0.030, MC statistics 0.021)
```

### `output/<src>_spectrum.png`

Left: the spectrum at the requested FCCD, deposited (grey), after the dead layer (blue) and after the
resolution (red), log scale, 0.2, 0.5 or 1 keV bins. Right: the observable against the FCCD with its
MC error; with a measurement, its band (red) and the FCCD read off (star).

![Am-241](output/am241_spectrum.png)

Am-241: one line, so nothing above the peak. At 1 mm, more than half the energy is lost to the dead
layer. The shelf below the peak is 59.5 keV gammas absorbed in the transition layer, keeping only part
of their energy; before the dead layer (grey) it is nearly empty, since a Compton scatter of 59.5 keV
leaves at most 11.3 keV. The spike below 2 keV is gammas absorbed in the dead layer. The peak falls
by a factor 23 from 0 to 3 mm.

![Ba-133](output/ba133_spectrum.png)

Ba-133: the 79.6+81.0 keV doublet loses far more to the dead layer than the 356.0 keV peak, which is
the whole method: R falls from 1.21 to 0.38 over 3 mm. The 30-36 keV Cs X-rays and the 437 keV sum
(81.0 + 356.0) are there because it decays; the X-rays only in grey, since 1 mm of dead layer stops
all but 0.1% of them.

![Co-60](output/co60_spectrum.png)

Co-60: the two peaks, their Compton edges, and the 2505.7 keV sum peak (both gammas in the crystal,
1.29 × 10⁻⁴ per decay), which only decays give; near 75 keV, lead X-rays from the castle. The dead
layer acts through the volume it removes: the 1332.5 keV peak falls about 5% per 0.25 mm.

![Th-228](output/th228_spectrum.png)

Th-228: the chain in equilibrium: Pb-212's 238.6, Tl-208's 583.2, 860.6 and 2614.5, Bi-212's 727.3
and 1620.5 keV; Pb and Bi X-rays at 73-88 keV and the 511 keV annihilation peak; the 2614.5 keV
Compton edge at 2382 keV, its single (2103.5) and double (1592.5) escape peaks, and sums above
2615 keV.

## Results

2 × 10⁶ to 5 × 10⁶ primaries per source, the public dummy `V99000A`, FCCD 1 mm, DLF 0.5 (the full
scans are in the printouts):

| | Am-241 | Ba-133 | Co-60 | Th-228 |
| --- | --- | --- | --- | --- |
| energy lost to the dead layer | 54.4% | 20.1% | 7.7% | 8.5% |
| observable | 59.5 keV: 7.69 × 10⁻⁶ per decay | R = 0.799 ± 0.007 | 1332.5 keV: 9.82 × 10⁻³ per decay | 2614.5 keV: 1.48 × 10⁻³ per decay |
| MC statistics on it | 1.0% | 0.8% | 0.7% | 1.2% |
| **a 1% error on it moves the FCCD by** | **0.010 mm** | **0.026 mm** | **0.051 mm** | **0.050 mm** |

The last row is each source's power. Am-241 is the sharpest dead-layer probe but needs the source
activity; Ba-133 needs none; Co-60 and Th-228 see the dead layer only through the volume it removes,
which is why Co-60 measures the active volume.

**Closure tests.**

- **The inversion:** R at 0.80 mm, fed back ± 0.01, returns 0.794 ± 0.037 mm.
- **The resolution:** Th-228's 583-2615 keV peaks return A = 0.867 ± 0.052 keV² (0.865 injected) and
  B = (2.255 ± 0.065) × 10⁻³ keV (2.250 × 10⁻³); each fitted width is within about 2σ of the injected
  one.

### Line mode against decays

The `am241` command with other cones, more gammas, or decays, each through `spectrum.C` (files of
0.3-2 GB, not kept): `THETA=3`; `NEV=44000000 SEED=5`; `THETA=180 NEV=50000000 SEED=6`; and
`PARTICLE=ion E=0 THETA=180 NEV=50000000 SEED=7`. The 59.5 keV peak:

| simulated as | primaries | per decay, FCCD 0 [10⁻⁵] | per decay, FCCD 1 mm [10⁻⁶] | FWHM at 1 mm [keV] |
| ------------ | --------- | ------------------------ | --------------------------- | ------------------ |
| the line, 3° cone | 4 × 10⁶ | 2.190 ± 0.004 | 7.624 ± 0.022 | 1.003 ± 0.003 |
| the line, 10° cone: `output/am241.root` | 4 × 10⁶ | 2.199 ± 0.012 | 7.688 ± 0.074 | 1.018 ± 0.009 |
| the line, 10° cone | 4.4 × 10⁷ | 2.225 ± 0.004 | 7.756 ± 0.022 | 1.005 ± 0.002 |
| the line, 4π | 5 × 10⁷ | 2.283 ± 0.041 | 7.67 ± 0.24 | 0.964 ± 0.031 |
| Am-241 decays | 5 × 10⁷ | 2.33 ± 0.07 | 8.66 ± 0.42 | 1.091 ± 0.064 |

Against the 4.4 × 10⁷ run: 3°, just wider than the channel's 2.24°, is 1.6 ± 0.2% low (1.7 ± 0.4% at
1 mm), missing the gammas that cross the collimator's copper edges. 4π (+2.6 ± 1.8%) and decays
(+4.8 ± 3.1%) agree with 10° within 1.5σ at FCCD 0, so `am241` runs at 10°; at 1 mm the decays are
2σ high, and their statistics cannot test the cone at the 1% level. The run in `output/` is
1.2 ± 0.6% low at FCCD 0 and 0.9 ± 1.0% at 1 mm, and its 1.018 keV width is high: both are
fluctuations, since with 11 times the gammas the width is 1.005 ± 0.002 keV (0.999 injected).

## On a measurement

1. Count the measured peaks in the windows `spectrum.C` prints: net = window − side bands.
2. Form the observable. Ba-133: R = net(79.6+81.0) / net(356.0), no activity needed. Otherwise net
   counts / (source activity × live time), per decay.
3. Read it back as an FCCD; the error combines the measurement and the MC statistics, not the
   limits below:

```bash
root -l -b -q 'ana/spectrum.C("output/ba133.root", 1.0, ‹R›, ‹error on R›)'
```

## Limits

- **The dead layer** is one FCCD with a linear transition (DLF 0.5), measured from the nearest
  surface. `dist_to_surf` stands in for the distance to the n+ contact, which differs near the p+
  contact and the groove, so the FCCD is an *effective* one.
- **The detector** is the public dummy `V99000A` in a dummy cryostat.
- **The resolution** is a placeholder, not a calibration.
- **Events** are one detector's steps within 10 µs: no threshold, pile-up or dead time.
- **Errors** are MC statistics only, one seed per run.

## Future work

| ‹placeholder› | adds | needs | plugs into |
| ------------- | ---- | ----- | ---------- |
| ‹field map› | A/E: how much of the Th-228 DEP, SEP and FEP survives | the impurity profile → field and weighting potential → drift-time map | a response step after remage (reboost's `drift_time`, `maximum_current`) |
| ‹pulse shapes› | waveforms: A/E with noise, rise times, surface-event cuts | the field map, the electronics response, measured noise | the same step |
| ‹real detector› | the crystal's real shape and cryostat | `legend-metadata` and the HADES cryostat records | `detector:` in `sim/<src>.yaml`, without `--public-geom` |
| ‹n+ contact› | a physical, not effective, FCCD | the distance to the n+ contact | the distance fed to `activeness()` |
| ‹transition layer› | a measured charge-collection profile | a measured or modelled profile | `activeness()` |
| ‹resolution› | spectra comparable to data | the detector's FWHM(E) | `cfg::resoA`, `cfg::resoB` |
| ‹activity and live time› | absolute efficiencies: Am-241 FCCD, Co-60 active volume | source certificates, DAQ live time | the `measured` argument |
| ‹active volume in cm³› | the active volume itself | the crystal volume and its surface map | from the FCCD Co-60 returns |
| ‹surface scans› | the FCCD across the surface | the campaign's source positions | `source_position` in `sim/<src>.yaml`, one run each |
| ‹DAQ› | threshold, pile-up, dead time | the DAQ settings | event building in `spectrum.C` |
| ‹systematics› | errors beyond MC statistics | runs with varied seed, source position, geometry | rerun, compare the FCCD |
