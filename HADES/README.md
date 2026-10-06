# HPGe characterization at HADES

One HPGe in a vacuum cryostat with a calibration source on top: the stand LEGEND characterizes its
detectors in, built by [legend-pygeom-hades](https://legend-pygeom-hades.readthedocs.io) and simulated
with remage. **The unknown is the detector, not the background**: each source is a known input, and
what the crystal makes of it measures its dead layer (FCCD), its active volume and its resolution.

| source | stand    | measures      | from                                              | simulated as                |
| ------ | -------- | ------------- | ------------------------------------------------- | --------------------------- |
| Am-241 | `am_HS1` | dead layer    | the 59.5 keV peak per decay (collimated)          | the 59.5 keV line, aimed    |
| Ba-133 | `ba_HS4` | dead layer    | R = 79.6+81.0 keV peaks / 356.0 keV peak: no activity needed | decays           |
| Co-60  | `co_HS5` | active volume | the 1332.5 keV peak per decay                     | decays                      |
| Th-228 | `th_HS2` | resolution, PSD peaks | lines from 239 to 2615 keV, DEP and SEP of 2614.5 keV | decays, the whole chain |

**Contents.** [Files](#files) · [Run](#run) · [How it works](#how-it-works) · [Results](#results) ·
[On a measurement](#on-a-measurement) · [Assumptions](#assumptions) · [Advanced: placeholders](#advanced-placeholders)

## Files

| file | |
| ---- | --- |
| `am241.yaml` `ba133.yaml` `co60.yaml` `th228.yaml` | the stand for each source: detector, `measurement:`, source position |
| `run.mac`    | remage: decay the source (or aim one gamma line) and record every step in the crystal |
| `spectrum.C` | steps → dead layer → resolution → peaks → FCCD scan → a measurement turned into an FCCD |
| `output/`    | everything generated: `<src>.gdml`, `<src>.root`, `<src>_spectrum.png` (git keeps the PNGs) |

No Python of our own: legend-pygeom-hades reads its YAML, the rest is remage and ROOT.

## Run

```bash
~/venvs/v/bin/pip install legend-pygeom-hades
```

remage v0.26.0 (Geant4 11.4.2) and ROOT 6.40 are the only other things needed. From `HADES/`:

```bash
mkdir -p output
~/venvs/v/bin/legend-pygeom-hades --public-geom --config am241.yaml output/am241.gdml
~/venvs/v/bin/legend-pygeom-hades --public-geom --config ba133.yaml output/ba133.gdml
~/venvs/v/bin/legend-pygeom-hades --public-geom --config co60.yaml output/co60.gdml
~/venvs/v/bin/legend-pygeom-hades --public-geom --config th228.yaml output/th228.gdml
remage -q --ignore-warnings -t 8 -w -o output/am241.root -s GDML=output/am241.gdml -s PARTICLE=gamma -s Z=95 -s A=241 -s E=59.5409 -s THETA=10 -s NEV=4000000 -s SEED=1 -- run.mac
remage -q --ignore-warnings -t 8 -w -o output/ba133.root -s GDML=output/ba133.gdml -s PARTICLE=ion -s Z=56 -s A=133 -s E=0 -s THETA=180 -s NEV=2000000 -s SEED=4 -- run.mac
remage -q --ignore-warnings -t 8 -w -o output/co60.root -s GDML=output/co60.gdml -s PARTICLE=ion -s Z=27 -s A=60 -s E=0 -s THETA=180 -s NEV=2000000 -s SEED=3 -- run.mac
remage -q --ignore-warnings -t 8 -w -o output/th228.root -s GDML=output/th228.gdml -s PARTICLE=ion -s Z=90 -s A=228 -s E=0 -s THETA=180 -s NEV=5000000 -s SEED=2 -- run.mac
root -l -b -q 'spectrum.C("output/am241.root")'
root -l -b -q 'spectrum.C("output/ba133.root")'
root -l -b -q 'spectrum.C("output/co60.root")'
root -l -b -q 'spectrum.C("output/th228.root")'
```

| source | primaries | remage, `-t 8` | file |
| ------ | --------- | -------------- | ---- |
| Am-241 | 4 × 10⁶ gammas, 10° cone | 31 s | 173 MB |
| Ba-133 | 2 × 10⁶ decays | 86 s | 155 MB |
| Co-60 | 2 × 10⁶ decays | 114 s | 161 MB |
| Th-228 | 5 × 10⁶ decays | 261 s | 281 MB |

Each GDML takes ~7 s and each `spectrum.C` ~5 s. More statistics: raise `NEV`, or add runs with
other seeds.

`--public-geom` builds the package's public dummy detector and cryostat, so no `legend-metadata`
checkout is needed. `spectrum.C` takes more arguments when you need them:

```
spectrum.C(file, fccd = 1.0, measured = -1, error = 0, dlf = 0.5, seed = 1)
  fccd       the FCCD [mm] the spectrum and the peak table are drawn at
  measured   a measured observable, with its error: printed back as an FCCD
  dlf        the fraction of the FCCD that is fully dead
  seed       the resolution smearing
```

## How it works

**`run.mac`** fires the source in one of two modes, chosen by its aliases:

| mode  | aliases                                    | what Geant4 does                                       | used for |
| ----- | ------------------------------------------ | ------------------------------------------------------ | -------- |
| decay | `PARTICLE=ion Z A E=0 THETA=180`           | the nuclide decays: every line, X-ray and coincidence sum | Ba-133, Co-60, Th-228 |
| line  | `PARTICLE=gamma E=<keV> THETA=<deg>`       | one gamma energy, emitted into a cone of `THETA` around the crystal | the collimated Am-241 |

- The source is confined to `Source_PV`, the disc that decays (`source_pv` is the holder around it).
- A daughter living longer than a year stops the decay (`DaughterNucleusMaxLifetime 1 y`): Np-237
  ends Am-241, while the Th-228 chain (members up to 3.6 d) runs to Pb-208, as in a source in equilibrium.
- remage stores every step in the crystal with its distance to the surface, and every primary in
  `stp/particles`: that is how `spectrum.C` knows the source, the mode and the cone.
- Line mode exists because the `am_HS1` channel (1.0 mm wide, 25.6 mm long) passes only 3.8 × 10⁻⁴ of
  the gammas: decays would waste 99.96% of the time. Aiming is checked against the real thing
  (5 × 10⁷ each): per emitted 59.5 keV gamma the full-energy peak is 6.22 × 10⁻⁵ in the 10° cone,
  6.25 × 10⁻⁵ ± 1.8% for gammas in 4π and 6.42 × 10⁻⁵ ± 3% for Am-241 decays. A 3° cone, just wider
  than the channel's 2.24°, misses 1.0 ± 0.2% that leaks through the collimator's edges; 10° does not.

**`spectrum.C`**, section by section:

```
1. configuration  resolution FWHM(E) = sqrt(A + B E), A 0.865 keV², B 0.00225 keV (a placeholder),
                  the 10 µs event window, the FCCD scan 0-3 mm, each source's peaks and observable
2. the source     from stp/particles. a nucleus: counts are per decay. a gamma: per decay = per gamma
                  x cone fraction (1 - cos THETA)/2 x the line's emission probability
3. events         steps sorted by event and time; a new detector event after 10 µs (splits the Th-228
                  chain, whose members decay days apart). deposited energy, and the active energy at
                  every FCCD at once: activeness 0 in the dead part (DLF x FCCD), linear to 1 at the FCCD
4. peaks          one Gaussian draw per event, shared by every FCCD; net counts = a ±3 sigma window
                  minus side bands of half its width on each side, as on a measured spectrum
5. response       at the requested FCCD: per peak the net counts, per decay, and the fitted FWHM against
                  the injected one; FWHM² = A + B E fitted across the peaks gives the resolution back
6. FCCD scan      the observable at every FCCD; a measured value ± error read back as FCCD ± error
7. draw           output/<src>_spectrum.png: the spectrum (deposited, active, smeared) | the scan
```

## Results

2 × 10⁶ to 5 × 10⁶ primaries per source, the public dummy `V99000A`, at FCCD 1 mm and DLF 0.5:

| | Am-241 | Ba-133 | Co-60 | Th-228 |
| --- | --- | --- | --- | --- |
| energy lost to the dead layer | 54.4% | 20.1% | 7.7% | 8.5% |
| observable | 59.5 keV: 7.69 × 10⁻⁶ per decay | R = 0.799 ± 0.007 | 1332.5 keV: 9.82 × 10⁻³ per decay | 2614.5 keV: 1.48 × 10⁻³ per decay |
| MC statistics on it | 1.0% | 0.8% | 0.7% | 1.2% |
| **a 1% error on it moves the FCCD by** | **0.010 mm** | **0.026 mm** | **0.051 mm** | **0.050 mm** |

The last row is each source's power: Am-241 at 59.5 keV is the sharpest dead-layer probe but needs
the source activity; Ba-133 needs no activity; Co-60 and Th-228 see the dead layer only through the
volume it removes, which is why Co-60 is the active-volume measurement.

### The FCCD scan

The observable against the dead layer: a measurement read off this table (`spectrum.C` interpolates
it) is the detector's FCCD.

| FCCD [mm] | Am-241: 59.5 keV per decay [10⁻⁶] | Ba-133: R | Co-60: 1332.5 keV per decay [10⁻³] | Th-228: 2614.5 keV per decay [10⁻³] |
| --- | --- | --- | --- | --- |
| 0.00 | 21.99 ± 0.12 | 1.213 ± 0.008 | 11.84 ± 0.08 | 1.787 ± 0.019 |
| 0.25 | 16.60 ± 0.11 | 1.082 ± 0.008 | 11.29 ± 0.08 | 1.697 ± 0.019 |
| 0.50 | 12.79 ± 0.10 | 0.975 ± 0.007 | 10.77 ± 0.07 | 1.621 ± 0.018 |
| 0.75 | 9.92 ± 0.08 | 0.882 ± 0.007 | 10.30 ± 0.07 | 1.550 ± 0.018 |
| **1.00** | 7.69 ± 0.07 | 0.799 ± 0.007 | 9.82 ± 0.07 | 1.478 ± 0.017 |
| 1.25 | 5.98 ± 0.07 | 0.727 ± 0.006 | 9.34 ± 0.07 | 1.402 ± 0.017 |
| 1.50 | 4.63 ± 0.06 | 0.658 ± 0.006 | 8.89 ± 0.07 | 1.335 ± 0.016 |
| 1.75 | 3.56 ± 0.05 | 0.599 ± 0.006 | 8.49 ± 0.07 | 1.271 ± 0.016 |
| 2.00 | 2.70 ± 0.04 | 0.544 ± 0.006 | 8.08 ± 0.06 | 1.204 ± 0.016 |
| 2.50 | 1.61 ± 0.04 | 0.451 ± 0.005 | 7.31 ± 0.06 | 1.086 ± 0.015 |
| 3.00 | 0.97 ± 0.03 | 0.375 ± 0.005 | 6.59 ± 0.06 | 0.977 ± 0.014 |

### Ba-133, as `spectrum.C` prints it

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

![Ba-133](output/ba133_spectrum.png)

Left: the spectrum as deposited (grey), after the 1 mm dead layer (blue) and after the resolution
(red); the 79.6+81.0 keV doublet loses far more to the dead layer than the 356.0 keV peak, which is the
whole method. Right: R against the FCCD.

### Am-241, Co-60, Th-228

![Am-241](output/am241_spectrum.png)

The aimed 59.5 keV line: more than half its energy lands in dead material at 1 mm, and the peak per
decay falls by a factor 23 from 0 to 3 mm. The fitted width, 1.018 ± 0.009 keV against 0.999
injected, is broadened by the low-energy shoulder of events that end in the transition layer.

![Co-60](output/co60_spectrum.png)

1173.2 keV: 1.105 × 10⁻² per decay, 1332.5 keV: 9.82 × 10⁻³, and the 2505.7 keV sum peak at
1.29 × 10⁻⁴: both gammas in the crystal at once, which only a decay simulation gives.

![Th-228](output/th228_spectrum.png)

The whole chain in equilibrium, per Th-228 decay: 2614.5 keV 1.48 × 10⁻³, its single escape (SEP,
2103.5 keV) 2.09 × 10⁻⁴, its double escape (DEP, 1592.5 keV) 9.5 × 10⁻⁵, and Bi-212's 1620.5 keV
9.9 × 10⁻⁵: the populations pulse-shape discrimination is calibrated on (‹field map›, below). Fitted
across 583-2615 keV, the peak widths give the resolution back: A 0.867 ± 0.052 keV² (0.865 injected),
B (2.255 ± 0.065) × 10⁻³ keV (2.250 × 10⁻³).

## On a measurement

1. Count the peaks of the measured spectrum with the windows `spectrum.C` prints (net = window −
   side bands).
2. Form the observable. Ba-133: R = net(79.6+81.0) / net(356.0), no activity needed. Am-241 and
   Co-60: net counts / (source activity × live time), per decay.
3. Read it back as an FCCD:

```bash
root -l -b -q 'spectrum.C("output/ba133.root", 1.0, ‹R›, ‹error on R›)'
```

Closure: R at FCCD 0.80 mm (0.867), fed back as a measurement ± 0.01, returns 0.794 ± 0.037 mm.

## Assumptions

- **The dead layer** is one FCCD with a linear transition (DLF 0.5), measured from the nearest
  surface. remage's `dist_to_surf` stands in for the distance to the n+ contact; they differ near
  the p+ contact and the passivated groove, so the FCCD here is an *effective* one.
- **The detector is the public dummy** `V99000A` in a dummy cryostat, not a real crystal.
- **The resolution is a placeholder**, not this detector's calibration.
- **Events** are one detector's steps within 10 µs; no trigger threshold, pile-up or dead time.
- **One seed per run**: the errors are MC statistics only.

## Advanced: placeholders

Each `‹placeholder›` is a piece this simulation does not do yet: what it adds, what it needs, and
where it plugs in.

| ‹placeholder› | adds | needs | plugs into |
| ------------- | ---- | ----- | ---------- |
| ‹field map› | A/E and pulse-shape discrimination: how much of the Th-228 DEP, SEP and FEP survives | the crystal's impurity profile → electric field and weighting potential (fieldgen/siggen or SolidStateDetectors.jl) → a drift-time map | a response step after remage: reboost's `drift_time` and `maximum_current` give A/E per event |
| ‹pulse shapes› | waveforms: A/E with noise, rise times, surface-event cuts | the field map, the electronics response, measured noise | the same response step |
| ‹real detector› | the crystal's real shape and cryostat | `legend-metadata` (diode) and the HADES cryostat records | `detector:` in the YAML; drop `--public-geom` |
| ‹n+ contact› | a physical, not effective, FCCD | the distance to the n+ contact, not to the nearest surface | the distance fed to `activeness()` |
| ‹transition layer› | a measured charge-collection profile instead of the linear one | a measured or modelled profile | `activeness()` |
| ‹resolution› | spectra comparable to data | the detector's calibration, FWHM(E) | `cfg::resoA`, `cfg::resoB` |
| ‹activity and live time› | absolute efficiencies: Am-241 FCCD, Co-60 active volume | source certificates, DAQ live time | the `measured` argument |
| ‹active volume in cm³› | the active volume itself, not the FCCD equivalent | the crystal volume and its surface map | from the FCCD Co-60 returns |
| ‹surface scans› | the FCCD across the surface: `am_HS1` at other r and φ, `th_HS2` lateral | the campaign's source positions | `source_position` in the YAML, one run per position |
| ‹DAQ› | trigger threshold, pile-up, dead time | the DAQ settings | event building in `spectrum.C` |
| ‹systematics› | errors beyond MC statistics | runs with varied seeds, source position, geometry | rerun, compare the FCCD |
