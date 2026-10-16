# LEGEND-1000 re-entrance tube

The re-entrance tube (RT) holds the underground argon (UGLAr) and the detector array. It is a 6 m × 2 m cylinder of
steel (SS) at the top, oxygen-free (OFHC) copper (Cu) in the middle and electroformed copper (EFCu) around the detectors.

**The question:** how much of the expensive EFCu can be cut, and how much cheap SS added, without changing the
tube's background index (BI)? "Without changing" means within 10%: that is the BI's own statistical uncertainty
(±10%), so a smaller change cannot be told apart from none.

**The answer** (depths below the tube's top; as built the SS ends at 2.06 m and the EFCu starts at 4.07 m):

|                                                                 | SS : Cu : EFCu     | BI after all cuts [cts/(keV·kg·yr)]          |
| --------------------------------------------------------------- | ------------------ | -------------------------------------------- |
| as built                                                        | 619 : 672 : 175 kg | 0.98 ± 0.10 × 10⁻⁸ (before cuts 4.60 × 10⁻⁶) |
| EFCu from 4.50 m (0.43 m shorter), SS to 2.28 m (0.22 m longer) | 685 : 635 : 137 kg | +10%: within its own ±10%                    |

[how radioactive is EFCu?](#how-radioactive-is-efcu) · [Method](#method) · [Results](#results) ·
[Future work](#future-work)

## how radioactive is EFCu?

What electroformed copper's own radioactivity looks like in the detectors.

`hpge_string_support_weldment_copper`: the 1008 (1.92 kg) EFCu detector holders, clamping the PEN plate under the detector and holding it to support rods, are simulated with EFCu's Th-232 and U-238 chains.

Two lines threaten `Qββ`: **Tl-208 at 2614 keV** and
**Bi-214 at 2448 keV**, whose gammas can leave 2039 keV in a detector.

Th-chain: 370000 chains, 5005703 germanium steps -> 434263 detector events

U-chain: 190000 chains, 1960954 germanium steps -> 176339 detector events

- **Simulated** (`sim/holders.mac`): 370 k Th-232 and 190 k U-238 chains, the ratio of EFCu's activities (0.37 : 0.19 µBq/kg), so the two spectra add as counted.
- **The spectrum** (`ana/holders.C`): raw energy each of 336 detectors gets per event.

```bash
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --holders output/holders_geometry.png
```

![the holders](output/holders_geometry.png)

```bash
remage -q --ignore-warnings -t 8 -w -o output/holders_th232.root -s GDML=geom/l1000geom/l1000.gdml -s Z=90 -s A=232 -s NEV=370000 -s SEED=21 -- sim/holders.mac
remage -q --ignore-warnings -t 8 -w -o output/holders_u238.root -s GDML=geom/l1000geom/l1000.gdml -s Z=92 -s A=238 -s NEV=190000 -s SEED=22 -- sim/holders.mac
root -l -b -q ana/holders.C
```

![the holders' spectrum](output/holders_spectrum.png)

`sim/holders.mac` writes `output/holders_th232.root` and `holders_u238.root`: `stp/germanium`, every germanium step
(`evtid`, `det_uid`, `edep_in_keV`, `time_in_ns`, ...), and `stp/vtx`, one row per chain; `ana/holders.C` draws them.

## Method

**Quick run**, from the NERSC decays already in `output/output_lh5`: the two response maps, the detector response,
the analysis.

```bash
rsync -a --partial nersc:/global/cfs/cdirs/m2676/users/neuberger/L1000_optical_muon_sims/omaps/v0.4.0/ular/merged/merged_optmap_20260225_063750.lh5 output/
rsync -a nersc:/global/cfs/cdirs/m2676/users/neuberger/L1000_optical_muon_sims/hpge_related/dtmaps/gen/V00000A.lh5 output/dtmap_V00000A.lh5
~/venvs/v/bin/python sim/hit.py output/output_lh5/*.lh5
root -l -b -q ana/background.C
```

`output/output_lh5` and `output/output_root` are links to `~/Downloads/output_lh5` and `~/Downloads/output_root`,
out of iCloud; git keeps only the PNGs in `output/`.

### 1. The geometry

The full LEGEND-1000 geometry, `geom/l1000geom/l1000.gdml` (legend-pygeom-l1000 0.6.0).

The tube is 6.246 m tall and 0.999 m in radius; the 336 detectors sit 4.8-5.9 m below its top, inside the
EFCu section. Walls: EFCu 1.5 mm (2.1 mm at the Cu seam), Cu and SS 6 mm; `ana/background.C` reads and checks it,
and its printout opens with it as TUBE: the outline, the seams at z = 920 and 2925 mm (4.07 m and 2.06 m below
the top, the depths a design moves), the wall thicknesses, the densities and the germanium mass:

```
TUBE  read from geom/l1000geom/l1000.gdml: z -1259.0 .. 4987.0 mm (6.246 m), 182 outline points, seams EFCu | 920.0 | Cu | 2925.0 | SS
    clean surface: on-axis points 2 (expect 2), duplicates 0, spikes 0, self-intersections 0 (expect 0)
    endcap: radius 999.0000 at z -946.6, barrel 999.0000: step +0.0000 mm
    shells against the barrel: Cu -0.0000, SS -0.0000 mm
    wall thickness [mm]:  EFCu 1.50 (z -947)  EFCu 1.50 (z -13)  EFCu 2.10 (z 910)  Cu 6.00 (z 1922)  SS 6.00 (z 3956)
    densities [kg/m^3], the GDML's: SS 8000, Cu 8960, EFCu 8960
    germanium, the GDML's: 336 detectors, 1023.4 kg: the BI is per kg of it
    tube: PASS
```

![the tube](output/tube.png)

The tube's three sections are the decay sources (SS : Cu : EFCu = 619.4 : 671.7 : 174.6 kg).

![the tube, the decay source](output/tube_geometry.png)

The underground argon around the detectors (18.83 m³, 26.17 t) is the veto, its light read by 12096 fibres and
252 SiPMs.

![the UGLAr, the veto](output/uglar_geometry.png)

Each detector is an inverted-coaxial point-contact (ICPC) germanium detector, 90.0 mm tall and 88.8 mm across, its
borehole 53 mm deep from the top: 3.046 kg, 1023.4 kg in all.

![the HPGe detectors](output/hpge_geometry.png)

`geom/l1000geom/l1000-autopeel.py` is the geometry viewer: `--holders`, `--tube`, `--uglar`, `--hpge <png>` draw one part
alone, offscreen, and print its numbers.

```bash
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --tube output/tube_geometry.png
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --uglar output/uglar_geometry.png
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --hpge output/hpge_geometry.png
```

### 2. Simulation

With `sim/run.mac` (remage), one nuclide decays uniformly through one section of the tube wall: **Tl-208** (²³²Th chain) or **Bi-214**
(²³⁸U chain), the only chain members whose gammas reach `Qββ`. It records the energy each decay leaves in each
detector and in the UGLAr. `sim/run.mac` takes `-s` aliases `GDML`, `Z` `A` (Tl-208: 81 208, Bi-214: 83 214),
`NEV`, `SEED`, and `VOLS`, a regex over the source volumes: `reentrance_tube_copper` (EFCu),
`reentrance_tube_layer_copper_ofhc` (Cu) or `reentrance_tube_layer_steel_316L` (SS).

Each run is 10 M decays (7.6 M in the EFCu runs) of one nuclide in one section, a 32-thread NERSC task of 16-60 min, and writes
`output/output_lh5/<nuclide>_<SEED>.lh5`: `stp/` (one table per germanium detector hit, and
`liquid_argon_underground`), `vtx` (one row per decay, with its position), `detector_origins`,
`number_of_simulated_events`. They are kept: every rerun starts here.

#### Why each section separately, and why these numbers

The three sections are simulated separately because their decays reach the detectors at very different rates, so each
can get the number of decays it needs; a run over the whole tube would put nearly 9 in 10 of its decays in the Cu and
SS, which the detectors barely see. §5 explains how the three are then combined into any design.

| section | runs | decays (Tl-208 / Bi-214) | real years they stand for | why that many                                                                                           |
| ------- | ---- | ------------------------ | ------------------------- | ------------------------------------------------------------------------------------------------------- |
| EFCu    | 6    | 23 M / 23 M              | 31 000 / 22 000           | next to the array, they give 94 of the BI's 113 events after all cuts: the BI's ±10%                    |
| Cu      | 366  | 460 M / 3.2 G            | 55 000 / 116 000          | they decide how far the SS can come down (§5): enough that one event, turned into SS, adds 3% of the BI |
| SS      | 36   | 50 M / 310 M             | 7 / 6                     | 2.7-4.8 m above the array: none of their decays puts an event in the ROI, at any cut level              |
| total   | 408  | 4.1 × 10⁹                |                           | 66 GB                                                                                                   |

"Real years" is how long the section, as built, takes to make that many decays of the nuclide (its activity × mass,
§4). In the Cu and SS, Bi-214 needs ~7 times the decays of Tl-208, because the SS holds ~7 times more Bi-214 than
Tl-208 (2500 vs 360 µBq/kg): for an event turned into SS to count the same in both chains, the two need the same
weight (§4), so 7 times the decays.

They ran on NERSC's Perlmutter: 122 task-hours, 17 node-hours of the m2676 allocation (16 min for a typical Cu run).
Back to back on this MacBook (M1 Pro, 8 threads: 22.5 min for the same run) they would take ~175 hours, about a week.
Another 55 node-hours went into runs that turned out not to change the answer, and were deleted.

To rerun them (`<SEEDS>`: EFCu 7001-7003 / 8001-8003, Cu 5001-5046 / 6001-6320, SS 3001-3005 / 4001-4031 for
Tl-208 / Bi-214; `<VOLS>` the section's volume, as in `sim/run.mac` above):

```bash
ssh nersc 'mkdir -p $SCRATCH/LEGEND1000-Simulation/sim $SCRATCH/LEGEND1000-Simulation/geom/l1000geom $SCRATCH/LEGEND1000-Simulation/output/output_lh5'
rsync -a geom/l1000geom/l1000.gdml nersc:'$SCRATCH/LEGEND1000-Simulation/geom/l1000geom/'
rsync -a sim/run.mac nersc:'$SCRATCH/LEGEND1000-Simulation/sim/'
ssh nersc 'shifterimg pull docker:legendexp/remage:v1.1.0'
ssh nersc 'cd $SCRATCH/LEGEND1000-Simulation && sbatch --account=m2676 --qos=shared --constraint=cpu --cpus-per-task=32 --time=01:30:00 --output=output/%x_%a.log --image=docker:legendexp/remage:v1.1.0 --array=<SEEDS> --job-name=rt --wrap "shifter --module=none remage -q --ignore-warnings -t 32 -w -m -s GDML=geom/l1000geom/l1000.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/output_lh5/<NUCLIDE>_\$SLURM_ARRAY_TASK_ID.lh5 -s Z=<Z> -s A=<A> -s VOLS=<VOLS> -- sim/run.mac"'
rsync -a --include='*_[0-9][0-9][0-9][0-9].lh5' --exclude='*' nersc:'$SCRATCH/LEGEND1000-Simulation/output/output_lh5/' output/output_lh5/
ssh nersc 'rm -r $SCRATCH/LEGEND1000-Simulation'
```

`SEED` starts remage's random numbers: the same seed repeats the same decays, so every run has its own, and it names
the run's file (`tl208_<SEED>.lh5`).

`sim/EdgarSim_analysis.py` and `sim/RalphSim_materialMix.ipynb` are reference only, for Edgar's older-geometry
simulation (`output/*_EFCu_RT_nolayer_optical_map_vtx_update.parquet`).

### 3. Detector response

`sim/hit.py` (reboost 1.4) gives, in the germanium, the energy through the dead
layer, the resolution and the A/E classifier; in the UGLAr, the photoelectrons the SiPMs see, from the
optical map. Its inputs: the response parameters of [legend1000-metadata](https://github.com/legend-exp/legend1000-metadata)
(`~/Documents/legend1000-metadata`, `simprod/config`, experiment `l1000dsg01`, as in
[legend-simflow](https://github.com/legend-exp/legend-simflow)), the UGLAr optical map v0.4.0
(`output/merged_optmap_20260225_063750.lh5`, its `/all` group, every SiPM summed) and the drift-time maps of `V00000A`
(`output/dtmap_V00000A.lh5`, 3500 V, ⟨100⟩ and ⟨110⟩ axes: how long each step's charge takes to reach the contact,
which shapes the current pulse and so A/E), applied with [reboost](https://reboost.readthedocs.io). Every detector is the dummy `V99999Z`; its drift times are those of `V00000A`, which has the same geometry (a stand-in: real detectors would each bring
their own dead layer, resolution and A/E).

|                     | value                                                                                                         | from                                                         |
| ------------------- | ------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| dead layer          | full-charge-collection depth (FCCD) 0.7 mm, outer half fully dead, linear in between                          | `V99999Z`; `tier/hit` `dead_layer_fraction` 0.5              |
| resolution          | FWHM = √(0.5 + 0.001·E) keV: 1.59 keV at `Qββ`                                                                | `pars/geds/eresmod`                                          |
| A/E                 | A: max of one template current pulse per step, at its drift time; σ_A/E = 0.01                                | `pars/geds/currmod`, `aoeresmod`; `output/dtmap_V00000A.lh5` |
| photoelectrons (PE) | LAr scintillation × the summed optical map; σ 0.3 PE, 16 ns merging, ≤ 100 PE per hit                         | `/all` of the optical map; `tier/opt`                        |
| M1, LAr veto        | one detector above 25 keV; rejects ≥ 4 PE summed                                                              | `tier/evt`                                                   |
| A/E cut             | classifier > −1.06, tuned to keep 90% of the Tl-208 double-escape peak (DEP); the production's −1.5 keeps 92% | LEGEND's convention                                          |

- **A/E** is the amplitude of the current-pulse over the energy. A gamma that scatters several times
  in the crystal spreads its charge out in time, giving a broad, low pulse: A/E << 1.
  0νββ puts its energy (two electrons) in one spot: A/E ≈ 1. The cut removes events with low
  A/E, and LEGEND places it so that 90% of the Tl-208 double-escape peak (1592.5 keV, single-site) survives.
  Here it is one value, tuned on simulation from a single current template, not a pulse-shape library.
  The printout's A/E CUT part: the DEP it is tuned on, and what the cut keeps of each Tl-208 peak (DEP; SEP, single
  escape; FEP, full energy) and of the ROI, next to `V00000A`'s measurement at HADES:

  ```
  A/E CUT  tuned on the Tl208 DEP (1592.5 keV) in M1 hits: 1025 in +-2.5 keV, 1073 in the sidebands 5-10 keV off it, 488 net
      classifier > -1.06 keeps 90% of it (+- 1% statistics); the production's > -1.50 kept 92%
      kept by it, net of the continuum:  DEP 1592.5 0.90 (488 net)  SEP 2103.5 0.02 (2810 net)  FEP 2614.5 0.04 (35578 net)  ROI 0.25 (42212)
      measured on V00000A at HADES:       DEP 1592.5 0.90            SEP 2103.5 0.04            FEP 2614.5 0.06            ROI 0.27
  ```

- **Optical map** (Moritz Neuberger's, on NERSC): for every 5 mm cell of the UGLAr, the probability that a
  scintillation photon made there reaches a SiPM. `sim/hit.py` turns each argon deposit into photons,
  and these into photoelectrons with that probability; the LAr veto counts them, all SiPMs summed (the
  production's second condition, 4 SiPMs with light, and random coincidences are left out). This is the v0.4.0
  map: the outer 68 mm of this UGLAr is read at its edge, and argon above z = 1.375 m gives no light, so until
  the 0.6.0 map arrives the BI is a 1.0-3.2 × 10⁻⁸ range (the BI AS BUILT lines read off the map alone give
  the top).

`sim/hit.py` takes remage outputs, plus `--optmap`, `--dtmap`, `--metadata`, `--gdml`, `--seed`, `--outdir`, `--force`
(defaults: the inputs above, seed 1, `output/output_root`). Its `--seed` (default 1) is not remage's `SEED`: with each file's name it starts the response's own random
numbers (energy smearing, photons, photoelectrons), so a file always gives the same hits. It writes
`output/output_root/<run>_hit.root`: TTrees `geds` (`evtid`, `det`, `edep`, `energy` [keV], `aoe_class`), `lar`
(`evtid`, `edep`, `pe`, `pe_map`), `vtx` (`evtid`, `xloc`, `yloc`, `zloc` [m]); ~50 GB. Only `ana/background.C`
reads them. They are kept: `sim/hit.py` skips a run whose hit file is newer than its inputs (`--force` rebuilds
it), and a rebuild gives the same hits, ~7 s a file.

With a new optical map, no new simulation: only the response and the analysis rerun.

```bash
~/venvs/v/bin/python sim/hit.py --optmap output/<new map>.lh5 output/output_lh5/*.lh5
root -l -b -q ana/background.C
```

### 4. Cuts and BI

`ana/background.C` is the analysis and its figures; it takes the tube runs (`file=isotope`, wildcards merge files)
and the GDML; the cuts, activities and the 10% line are in `cfg` (`m1_keV`, `lar_pe`, `depKeep`, `chain`, `simTol`).
Its printout, in order: TUBE (§1), A/E CUT (§3), INPUTS and, per chain, DECAYS (§4) and DEPTH (§5), BI AS BUILT and
DESIGNS (Results); each part is shown in full where it is used.

**Cuts**, in LEGEND's order: **M1** (exactly one detector above 25 keV), the **LAr veto** (rejects ≥ 4
photoelectrons), **A/E** (single-site-like, at the value that keeps 90% of the Tl-208 double-escape peak).

**BI** is counts in the **ROI**, MAJORANA's background region ([Haufe et al.](https://arxiv.org/abs/2209.10592),
[MAJORANA's final result](https://arxiv.org/abs/2207.07638)): 1950-2350 keV, less ±5 keV around each of 2039.0
(`Qββ`), 2103.5 (Tl-208 single escape), 2118.5 and 2204.1 keV (Bi-214), so **360 keV** of continuum, per keV, per kg
of germanium (1023.4 kg), per year.

**ROI vs `Qββ` ± 2σ.** The BI is a rate per keV, so it is counted over the wide 360 keV ROI to collect enough
events; the continuum there is flat, so its average holds at `Qββ`. The green band on the spectrum, `Qββ` ± 2σ
(2037.6-2040.4 keV, ~95% of a 0νββ peak at this resolution), is where a signal would sit: the background expected
under it is the BI × 2.7 keV × the exposure.

**Normalisation** turns simulated decays into real ones. Per section, from the GDML and the activities: M = ρV,
A = M × a (a: the chain's µBq/kg × the nuclide's share of it), N_real = A × 1 yr. The simulation put N_sim decays
there, so each stands for N_real / N_sim real decays a year (its weight). An ROI event adds that weight /
(1023.4 kg × 360 keV) to the BI. §5 uses these weights to try any design without a new simulation (bulk activities only,
chains in equilibrium).

**Activities µBq/kg** (`cfg::chain` in `ana/background.C`), from Ralph Massarczyk's MaterialMix talk (23 Jun 2026).
Bi-214 is 100% of the ²³⁸U chain; Tl-208 35.94% of the ²³²Th chain ([Nuclear Data Sheets](https://doi.org/10.1016/j.nds.2007.07.001)).

|      | ²³²Th chain | ²³⁸U chain | source quoted there                                           |
| ---- | ----------- | ---------- | ------------------------------------------------------------- |
| SS   | 1000        | 2500       | Bernhard                                                      |
| Cu   | 1.1         | 1.3        | [MAJORANA assay](https://doi.org/10.1103/PhysRevC.110.055804) |
| EFCu | 0.37        | 0.19       | M. Green                                                      |

`ana/background.C` prints the ROI, cuts and activities it uses as INPUTS:

```
INPUTS  ROI: 1950-2350 keV minus 10 keV around 2039, 2103.5, 2118.5, 2204.1 keV (MAJORANA's background region): 360 keV
        cuts: M1 (exactly one detector above 25 keV), LAr veto (the SiPMs see < 4 photoelectrons), A/E classifier > -1.06 (90% DEP, A/E CUT)
        activity: chain [uBq/kg] 232Th / 238U: SS 1000 / 2500, Cu 1.1 / 1.3, EFCu 0.37 / 0.19; Tl208 is 35.94% of 232Th, Bi214 100% of 238U
```

The printout's DECAYS part, per chain: where the decays fell, the germanium hits, and the normalisation (`simulated`
is N_sim, `per m^3` how densely the runs sampled the section, `weight` is N_real / N_sim; the sections are sampled
very unevenly on purpose, §2, so each uses its own density). The hit counts are cumulative: each row adds its cut to the ones above. They
are counted twice, over every energy and inside the ROI; the bottom-right number, in the ROI past all three cuts, is
what the BI is built from.

```
DECAYS  Tl208  (output/output_root/tl208*_hit.root, 532811152 decays, 54 files)
    where they decayed: EFCu 4.2%, Cu 86.4%, SS 9.4%; in the argon 0, outside the tube 0
    germanium hits, per cut level (each row adds its cut to the ones above), at all energies and inside the 360 keV ROI:
    cut level          all energies         in ROI
    no cut                  1507568          49788
    M1                      1288342          42212
    M1 + LAr                  32094            633
    M1 + LAr + A/E             4302             91
    LAr veto keeps 1 in 40 M1 hits; read off the optical map alone, without the argon beyond its edge, 1 in 14
    the normalisation, per section: volume, mass, activity, simulated decays (and per m^3), weight = real decays a year per simulated one
    section   V [m^3]   M [kg]     A [Bq]  simulated     per m^3   weight
    EFCu      0.01948    174.6  2.322e-05   22628745  1161362547 3.24e-05
    Cu        0.07497    671.7  2.656e-04  460181326  6138059588 1.82e-05
    SS        0.07742    619.4  2.226e-01   50001081   645805476     0.14
    the runs sampled the EFCu 0.342 x as densely as the shells (1 = evenly): each section's weight uses its own density

DECAYS  Bi214  (output/output_root/bi214*_hit.root, 3532803014 decays, 354 files)
    where they decayed: EFCu 0.6%, Cu 90.6%, SS 8.8%; in the argon 0, outside the tube 0
    germanium hits, per cut level (each row adds its cut to the ones above), at all energies and inside the 360 keV ROI:
    cut level          all energies         in ROI
    no cut                   838014           2684
    M1                       753211           2506
    M1 + LAr                  34223            601
    M1 + LAr + A/E             4567             22
    LAr veto keeps 1 in 22 M1 hits; read off the optical map alone, without the argon beyond its edge, 1 in 10
    the normalisation, per section: volume, mass, activity, simulated decays (and per m^3), weight = real decays a year per simulated one
    section   V [m^3]   M [kg]     A [Bq]  simulated     per m^3   weight
    EFCu      0.01948    174.6  3.317e-05   22619347  1160880219 4.63e-05
    Cu        0.07497    671.7  8.733e-04 3200182724 42685157225 8.61e-06
    SS        0.07742    619.4  1.548e+00  310000943  4003919566    0.158
    the runs sampled the EFCu 0.050 x as densely as the shells (1 = evenly): each section's weight uses its own density
```

### 5. Finding the best design: reweighting

A design moves the two seams: the SS ends at depth L1 (2.06 m as built) and the EFCu starts at depth L2 (4.07 m),
with Cu between. The wall keeps its shape and thickness as built; only the material of each part of it changes (so
EFCu turned into Cu at 4.07-4.50 m is the EFCu's 1.5 mm wall, made of Cu).

**Reweighting, not resimulating.** What a decay does in the detectors depends on where it happens, and very little on
what the wall there is made of: its gammas cross the same geometry either way, and SS and Cu stop them about equally
(8.0 vs 9.0 g/cm³). So each simulated decay can stand for any material a design puts at its depth; only its weight
(§4) changes: the real decays a year per m³ of that material (density × activity), over the simulated decays per m³
of the run it came from. A design's BI is the sum, over the events, of the weight each gets from the material the
design puts at its depth. Any L1 and L2 can be tried without a new simulation.

The decays at a depth always come from the section built there: the SS runs cover 0-2.06 m, the Cu runs 2.06-4.07 m,
the EFCu runs 4.07-6.25 m. So the two ways to save money reweight different runs:

| design change                     | decays reweighted                                  | weight × (Tl-208 / Bi-214) | events there                     |
| --------------------------------- | -------------------------------------------------- | -------------------------- | -------------------------------- |
| EFCu starts lower (L2 deeper)     | the EFCu runs', from 4.07 m to L2: EFCu becomes Cu | 3.0 / 6.8                  | thousands: the change is precise |
| SS comes further down (L1 deeper) | the Cu runs', from 2.06 m to L1: Cu becomes SS     | 810 / 1700                 | a handful: each one counts       |

**Why the Cu runs decide the SS.** The SS added in a design would sit where the Cu is now, 2-3 m above the array,
and only the Cu runs have decays there. Those decays rarely reach the detectors: one in 10⁷-10⁸ Tl-208 decays and one
in ~10⁹ Bi-214 decays gives an ROI event after M1 (DEPTH, below). Turned into SS, each such event counts ~1000 times
more. With 460 M / 3.2 G Cu decays, one of them adds 3.2% of the BI, so the 10% line is about three events: a
measurement, not one chance event. With only a tenth of them, one event would add ~30%; none happened to fall
between 2.06 and 2.66 m, and the SS looked free down to 2.66 m: the wrong answer.

**Cut survival.** Too few events survive all cuts away from the detectors to count the change on them, so it is
counted on the ROI events that pass M1, 100-500 times more, times the **cut survival**: the fraction of them that also
pass the LAr veto and A/E. It is measured where the change is:

- near the detectors (below 4.37 m), per chain: 1.9 × 10⁻³ (Tl-208), 8.6 × 10⁻³ (Bi-214);
- higher up, both chains together: 20 of 2559 = 7.8 × 10⁻³ (Tl-208 1 + 13 + 322 + 1813, Bi-214 2 + 64 + 344 = 2559
  events after M1, 20 after all cuts). Tl-208's events survive ~4 × more often there, partly because the optical
  map gives no light above z = 1.375 m.

Both come from the printout's DEPTH part, per 0.62 m slab below the tube's top: near, the rows below 4.37 m summed;
higher up, the rows above, both chains together:

```
DEPTH  Tl208  by depth below the tube's top: decays, germanium hits per decay, ROI events after M1 and after all cuts,
       and the cut survival = ROI events after all cuts / after M1
    depth [m]            decays   hits/decay      ROI: M1   all cuts   survival
     0.00..0.62       15146809     0.00e+00            0          0      -
     0.62..1.25       15140844     0.00e+00            0          0      -
     1.25..1.87       15148117     0.00e+00            0          0      -
     1.87..2.50      105323181     2.09e-07            1          0 0.0000
     2.50..3.12      144185058     3.93e-06           13          0 0.0000
     3.12..3.75      144192734     8.79e-05          322          6 0.0186
     3.75..4.37       74471336     8.67e-04         1813         10 0.0055
     4.37..5.00        6825938     4.53e-02         8796         18 0.0020
     5.00..5.62        6823148     1.07e-01        20243         42 0.0021
     5.62..6.25        5553987     6.98e-02        11024         15 0.0014   <- nearest the detectors

DEPTH  Bi214  by depth below the tube's top: decays, germanium hits per decay, ROI events after M1 and after all cuts,
       and the cut survival = ROI events after all cuts / after M1
    depth [m]            decays   hits/decay      ROI: M1   all cuts   survival
     0.00..0.62       93892792     0.00e+00            0          0      -
     0.62..1.25       93912243     0.00e+00            0          0      -
     1.25..1.87       93910381     0.00e+00            0          0      -
     1.87..2.50      729140823     8.23e-09            0          0      -
     2.50..3.12     1003152378     4.19e-07            2          0 0.0000
     3.12..3.75     1003111516     1.44e-05           64          1 0.0156
     3.75..4.37      496487697     1.75e-04          344          3 0.0087
     4.37..5.00        6819575     2.28e-02          445          4 0.0090
     5.00..5.62        6819893     5.59e-02         1077          9 0.0084
     5.62..6.25        5555716     3.59e-02          574          5 0.0087   <- nearest the detectors
```

**The scan.** For every L2, L1 goes down 1 cm at a time while the BI stays within 10% of the as-built's. The answer is
the least EFCu, then, for it, the most SS.

## Results

4066 M tube decays (Tl-208 / Bi-214): EFCu 23 M / 23 M, Cu 460 M / 3.2 G, SS 50 M / 310 M. After all cuts,
113 ROI events: 94 from the EFCu (76 Tl-208, 18 Bi-214), 19 from the Cu, none from the SS.

**The BI as built** (the LEGEND-1000 goal is 10⁻⁵):

| no cut      | M1                          | M1 + LAr                                     | all cuts (M1 + LAr + A/E)                           |
| ----------- | --------------------------- | -------------------------------------------- | --------------------------------------------------- |
| 4.60 × 10⁻⁶ | 3.92 × 10⁻⁶: M1 removes 15% | 1.11 ± 0.03 × 10⁻⁷: the LAr veto removes 97% | **0.98 ± 0.10 × 10⁻⁸**: A/E removes 91% of the rest |

The printout's BI AS BUILT part, per chain and section; `ROI hits` counts the ROI events at each cut level (no cut /
M1 / + LAr / + A/E):

```
BI AS BUILT  [cts/(keV kg yr)], +- simulation statistics; none: no ROI event in all its simulated decays
    chain  sect          ROI hits   no cut                 M1                     M1 + LAr               M1 + LAr + A/E
    Tl208  SS             0/0/0/0   none                   none                   none                   none
           Cu     1937/1751/176/15   9.57e-08 +- 2.2e-09    8.66e-08 +- 2.1e-09    8.70e-09 +- 6.6e-10    7.41e-10 +- 1.9e-10
           EFCu   47851/40461/457/76   4.20e-06 +- 1.9e-08    3.56e-06 +- 1.8e-08    4.02e-08 +- 1.9e-09    6.68e-09 +- 7.7e-10
           tube                     4.30e-06 +- 1.9e-08    3.64e-06 +- 1.8e-08    4.89e-08 +- 2.0e-09    7.42e-09 +- 7.9e-10
    Bi214  SS             0/0/0/0   none                   none                   none                   none
           Cu       414/396/131/4   9.68e-09 +- 4.8e-10    9.26e-09 +- 4.7e-10    3.06e-09 +- 2.7e-10    9.35e-11 +- 4.7e-11
           EFCu   2270/2110/470/18   2.85e-07 +- 6.0e-09    2.65e-07 +- 5.8e-09    5.90e-08 +- 2.7e-09    2.26e-09 +- 5.3e-10
           tube                     2.95e-07 +- 6.0e-09    2.74e-07 +- 5.8e-09    6.21e-08 +- 2.7e-09    2.35e-09 +- 5.3e-10
    ALL CHAINS
      no cut           4.60e-06 +- 2.0e-08 (0.46 x goal)
      M1               3.92e-06 +- 1.9e-08 (0.39 x goal)
      M1 + LAr         1.11e-07 +- 3.4e-09 (0.011 x goal)
      M1 + LAr + A/E   9.77e-09 +- 9.5e-10 (0.00098 x goal)
    ...with the LAr veto read off the optical map alone, no light from the argon beyond its edge:
      M1 + LAr         2.75e-07 +- 5.1e-09 (0.027 x goal)
      M1 + LAr + A/E   3.18e-08 +- 1.7e-09 (0.0032 x goal)
    ROI hits: per cut level, as in the columns. sections with no ROI hit are left out of the totals
```

![the spectrum](output/slide_spectrum.png)

The tube's spectrum as built, through each cut: the ROI shaded light green, `Qββ` ± 2σ dark green, the BI after all
cuts as the black line across the ROI, the 10⁻⁵ goal dashed.

**The designs**, for each depth L2 where the EFCu starts, the deepest L1 where the SS can end:

| EFCu from (L2)    | EFCu       | SS to (L1) | SS : Cu : EFCu         | BI change                      |
| ----------------- | ---------- | ---------- | ---------------------- | ------------------------------ |
| 4.07 m, as built  | 175 kg     | 2.69 m     | 809 : 460 : 175 kg     | +6%                            |
| 4.20 m            | 163 kg     | 2.69 m     | 809 : 472 : 163 kg     | +8%                            |
| 4.30 m            | 154 kg     | 2.66 m     | 800 : 490 : 154 kg     | +7%                            |
| 4.40 m            | 146 kg     | 2.66 m     | 800 : 499 : 146 kg     | +10%                           |
| **4.50 m**        | **137 kg** | **2.28 m** | **685 : 635 : 137 kg** | **+10%**                       |
| 4.60 m and deeper | ≤ 129 kg   | none       |                        | more than 10% at any SS length |

```
DESIGNS  against the tube as built (BI 9.77e-09 after all cuts): SS down to L1, Cu down to L2, EFCu below. a design passes when its BI
    stays within 10% of the as-built's, the BI's own statistical uncertainty. L1 from 2.06 m down, 1 cm at a time
    a design's BI change is counted on ROI events after M1 x the cut survival (see DEPTH): below 4.37 m Tl208 1.87e-03 Bi214 8.59e-03; above it, both chains, 20 of 2559 = 7.82e-03
    L2 [m]   EFCu | deepest L1 SS : Cu : EFCu                     BI change
    4.07   175 kg |     2.69 m  809 : 460 : 175 kg = 56 : 32 : 12 %       +6%   <- EFCu as built
    4.10   171 kg |     2.69 m  809 : 463 : 171 kg = 56 : 32 : 12 %       +7%
    4.20   163 kg |     2.69 m  809 : 472 : 163 kg = 56 : 33 : 11 %       +8%
    4.30   154 kg |     2.66 m  800 : 490 : 154 kg = 55 : 34 : 11 %       +7%
    4.40   146 kg |     2.66 m  800 : 499 : 146 kg = 55 : 35 : 10 %      +10%
    4.50   137 kg |     2.28 m  685 : 635 : 137 kg = 47 : 44 :  9 %      +10%
    4.60   129 kg | -          the BI rises more than 10% at any L1
    4.70   121 kg | -          the BI rises more than 10% at any L1
    4.80   112 kg | -          the BI rises more than 10% at any L1
    4.90   104 kg | -          the BI rises more than 10% at any L1
    5.00    95 kg | -          the BI rises more than 10% at any L1
    5.25    74 kg | -          the BI rises more than 10% at any L1
    5.50    53 kg | -          the BI rises more than 10% at any L1
    5.75    32 kg | -          the BI rises more than 10% at any L1
    6.00    11 kg | -          the BI rises more than 10% at any L1
    least EFCu, then most SS: L2 4.50 m, L1 2.28 m:  685 : 635 : 137 kg = 47 : 44 :  9 %, BI +10%
    as built:  619 : 672 : 175 kg = 42 : 46 : 12 %
```

- **EFCu** can start 0.43 m lower, at 4.50 m: 137 kg instead of 175 (−22%). Deeper, the Cu that would replace it
  sits next to the detectors at 3-7 × EFCu's activity, and the BI rises more than 10% whatever the SS.
- **SS** can then come down to 2.28 m: 685 kg instead of 619 (+66 kg). With the EFCu as built it could reach 2.69 m.
- **The SS as built** never reaches the detectors: no ROI event at any cut level in its 360 M decays.

![the designs](output/slide_designs.png)

The designs near the tube as built, coloured by their BI over the as-built's: the designs left of and below the bold
+10% line keep the BI. The star is the tube as built, the square the answer.

`/Applications/LEGEND.app`'s rt is the same scan (`output/designs.json`) on the GDML's tube in 3D: two sliders for
L1 and L2 recolour the tube by its BI and give the masses; `app/make_app.py` builds the app.

## Future work

**The Hall C geometry.** Matthew Busch's string and lock update (29 Sep 2026) proposes 12-detector strings, for
running beyond 1 t, at the same central depth in the UGLAr:

- **Array**: 42 strings × 8 detectors × 3 kg = 1008 kg becomes 48 × 12 × 3 kg = 1727 kg. Each string's germanium
  grows from 1.08 to 1.64 m and its fibres from 1.36 to 1.92 m, reaching 0.28 m higher and 0.28 m deeper.
- **Tube**: 0.2 m longer, the extra height in the Cu section, and a deeper bottom head (467 → 575 mm deep): more UGLAr.
- **Lab**: the cleanroom floor, cryostat and water tank drop 0.2 m; the lock pipe grows from 7.54 to 7.81 m. It
  still waits on the crane and cleanroom details.

He makes the change conditional on the background simulations, which is what this study can answer. The array
reaches 0.28 m higher, so the EFCu likely has to start higher too; the design scan says how much. It needs
the Hall C geometry from legend-pygeom-l1000 (`units: n: 12` in its `configs/string.yaml`, 48 strings, the longer
tube and head), then the same NERSC runs and analysis on it.
