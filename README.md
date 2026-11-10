# LEGEND 1000 Background Index

The re-entrance tube (RT) holds the underground argon (UGLAr) and the detector array (`LEGEND1000pCDR, LEGEND.pdf`,
p. 8). It is a 6 m × 2 m cylinder: steel (SS) at the top, oxygen-free high-conductivity (OFHC) copper (Cu) in the
middle, electroformed copper (EFCu) around the detectors (`ReentrantTubeFabrication, Lee.pdf`, p. 9). EFCu is 3-7 ×
less radioactive than Cu (the activities in step 4).

Sources are cited in place. Files named `….pdf` are in iCloud Drive, `0νββ/Papers` or `0νββ/Presentations`; "p." is
the PDF's page.

**The question:** how much of the EFCu can be cut, and how much SS added, without changing the tube's background
index (BI)? "Without changing" means within 10%:

1. Any material swap changes the BI a little, so "no change" can't mean exactly zero.
2. The as-built BI rests on 128 simulated counts after all cuts.
3. Simulated counts follow Poisson statistics (`ralph_GeneralCall_10_12_2026.pdf`, p. 20): a count of N scatters by
   √N, so √128 ≈ 11, ±9%.
4. So the BI itself is only known to ±9% (0.94 ± 0.09 × 10⁻⁸).
5. A change smaller than that can't be told apart from the simulation's own scatter.
6. So "no change" means within ~10%: a choice made for this study, not a LEGEND requirement, and the simulation was
   sized to resolve it (step 2, Simulation).

**The answer** (depths below the tube's top; as built the SS ends at 2.06 m and the EFCu starts at 4.07 m):

|                                                                 | SS : Cu : EFCu     | BI after all cuts [counts/(keV·kg·yr)]       |
| --------------------------------------------------------------- | ------------------ | -------------------------------------------- |
| as built                                                        | 619 : 672 : 175 kg | 0.94 ± 0.09 × 10⁻⁸ (before cuts 4.29 × 10⁻⁶) |
| EFCu from 4.40 m (0.33 m shorter), SS to 2.69 m (0.63 m longer) | 809 : 489 : 146 kg | +8%: within 10%                              |

[Getting started](#getting-started) · [How radioactive is EFCu?](#how-radioactive-is-efcu) · [Method](#method) ·
[Results](#results) · [LEGEND 1700?](#legend-1700)

## Getting started

- **Geometry:** legend-pygeom-l1000 0.6.0 writes the LEGEND-1000 geometry as GDML; remage reads it and registers the
  336 germanium detectors from the file itself.
- **Simulation:** remage 1.1.0 (`remage paper, LEGEND.pdf`) with one macro, `sim/sim.mac`: pixi on the Mac, the
  collaboration's container on NERSC.
- **Python:** `~/venvs/v` (Python 3.14: reboost 1.4.0, legend-pygeom-l1000 0.6.0, pylegendmeta, uproot, vtk).
- **Analysis:** ROOT 6.40.
- **Building from source** (Geant4, BxDecay0, remage; Apple Silicon and Linux; each stage tests itself):
  [Physics-Simulation-Stack](https://github.com/h0t5tuff/Physics-Simulation-Stack).
- **Not used:** the modified geometry `KSendcap_l1kGeometry.gdml`. Its tube is shifted 68 mm in radius and along the
  axis, and its endcap is 1 mm wider than the barrel where they meet, so it fails the tube check (step 1).

Three words, used throughout:

- **decay:** one simulated nuclear decay of Tl-208 or Bi-214 in the tube (a whole Th-232 or U-238 chain, for the holders);
- **hit:** the energy one decay leaves in one detector; one decay can leave hits in several detectors;
- **count:** a decay whose hit lands in the ROI (the 100 keV window around `Qββ` = 2039 keV) after the cuts. The
  simulation gives _simulated counts_; each stands for w real decays a year, and Σw is the _counts a year_.

## How radioactive is EFCu?

EFCu's own radioactivity, as the detectors see it:

- **Source:** the 1008 EFCu detector holders (`hpge_string_support_weldment_copper`, 1.92 kg in all; count and mass
  from the GDML), pictured below.
- **Activity:** 0.37 ± 0.22 µBq/kg (²³²Th chain) and 0.19 ± 0.10 µBq/kg (²³⁸U chain): the LEGEND-1000 CDR's Table
  XVII, shown in `CDRBackgroundModel, Green.pdf`, p. 7.
- **Decays:** 370 k Th-232 chains and 190 k U-238 chains (`sim/sim.mac`; `AMIN` `ZMIN` are the chain's lowest mass
  and atomic numbers). At these activities both are the holders' decays over the same ~16 500 years (1.92 kg ×
  activity × 1 yr = 22.4 and 11.5 chains a year), so the two spectra add as they are.
- **Spectrum:** `holders()` in `ana/ana.C` sums each detector's energy within 10 µs, with no detector response:
  435 530 hits from Th-232, 176 339 from U-238.
- **Lines:** two gamma lines are labelled, **Tl-208 at 2614 keV** and **Bi-214 at 2448 keV**
  (`Legend_Note_24_007.pdf`, p. 31). Both are above `Qββ` = 2039 keV, the energy of neutrinoless double-beta decay
  (0νββ) in ⁷⁶Ge (`LEGEND1000pCDR, LEGEND.pdf`, p. 7), which makes these two decays "of primary concern" (p. 29).

![the holders](output/holders_geometry.png)

![the holders' spectrum](output/holders_spectrum.png)

```bash
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --holders output/holders_geometry.png
remage -q --ignore-warnings -t 8 -w -o output/holders_th232.root -s GDML=geom/l1000geom/l1000.gdml -s VOLS='hpge_string_support_weldment_copper_.*' -s Z=90 -s A=232 -s ZMIN=81 -s AMIN=208 -s NEV=370000 -s SEED=21 -- sim/sim.mac
remage -q --ignore-warnings -t 8 -w -o output/holders_u238.root -s GDML=geom/l1000geom/l1000.gdml -s VOLS='hpge_string_support_weldment_copper_.*' -s Z=92 -s A=238 -s ZMIN=81 -s AMIN=206 -s NEV=190000 -s SEED=22 -- sim/sim.mac
root -l -b -q -e '.L ana/ana.C' -e 'holders()'
```

## Method

Five steps. Only the PNGs of `output/` are in git; the commands in each step remake the rest.

1. **Geometry:** build `geom/l1000geom/l1000.gdml` with legend-pygeom-l1000 and check it with
   `geom/l1000geom/check.mac` and `l1000-autopeel.py --table`, which writes the tube's numbers to `output/tube.txt`.
2. **Simulation:** `sim/sim.mac` (remage) decays Tl-208 and Bi-214 in each section of the tube on NERSC, into `output/output_lh5` (76 GB).
3. **Detector response:** `sim/hit.py` turns those decays into what the detectors and the LAr veto see, in
   `output/output_root`.
4. **Cuts and BI:** `ana/ana.C` applies the cuts and computes the BI.
5. **Reweighting:** `ana/ana.C` also tries every design, without a new simulation.

### 1. The geometry

Every number here is read from `geom/l1000geom/l1000.gdml` (legend-pygeom-l1000 0.6.0), the full LEGEND-1000 geometry:

- **Tube:** 6.246 m tall, 0.999 m in radius, its top at z = 4987 mm.
- **Seams,** the two a design moves: z = 920 mm (EFCu to Cu, 4.07 m deep) and z = 2925 mm (Cu to SS, 2.06 m deep).
- **Wall:** 1.5 mm of EFCu through the endcap and lower wall, thickening to 3 mm over the last 2 cm below the Cu
  seam; 6 mm through the Cu and SS, which are layers placed inside the EFCu tube volume.
- **Detectors:** 336 inverted-coaxial point-contact (ICPC) germanium detectors, 4.8-5.9 m below the top, inside the
  EFCu section. Each is 90.0 mm tall, 88.8 mm across, with a 53 mm borehole: 3.046 kg, 1023.4 kg in all.
- **Veto:** the UGLAr around the detectors (18.83 m³, 26.17 tonnes), its light read by 12096 fibres and 252 silicon
  photomultipliers (SiPMs).

Two checks run before any physics:

1. `check.mac` builds the geometry and scans it for overlapping volumes.
2. `l1000-autopeel.py --table` checks that the bottom head and the Cu and SS layers are flush with the barrel, then
   writes `output/tube.txt`: the seams, the densities, the germanium mass and the wall's cross-section every 0.5 mm.

Every later step reads the tube from `output/tube.txt` (the analysis, the pictures, the app), so a new design is just
a new GDML file. `l1000-autopeel.py` also draws each part alone, offscreen, and prints its numbers.

![the tube, the decay source](output/tube_geometry.png)

![the UGLAr, the veto](output/uglar_geometry.png)

![the HPGe detectors](output/hpge_geometry.png)

```bash
export LEGEND1000_METADATA=~/Documents/legend-exp/legend1000-metadata
~/venvs/v/bin/legend-pygeom-l1000 -v geom/l1000geom/l1000.gdml
remage --ignore-warnings -s GDML=geom/l1000geom/l1000.gdml -s SKIP=false -s NPOINTS=1000 -- geom/l1000geom/check.mac 2>&1 | grep "Overlap is detected"
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --table output/tube.txt
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --tube output/tube_geometry.png
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --uglar output/uglar_geometry.png
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --hpge output/hpge_geometry.png
```

**The lid**, the one part of the GDML this study leaves out:

1. In legend-pygeom-l1000 0.6.0 the argon closes 6 mm below the tube's top on the axis but at the top at the wall.
2. That leaves a 55.51 kg EFCu cone closing the tube's top (z = 4981-4987 mm): the lid, an artifact.
3. Here `sim/sim.mac` keeps it out of the sampling by its position.
4. [legend-exp/legend-pygeom-l1000#142](https://github.com/legend-exp/legend-pygeom-l1000/pull/142) fixes it: the
   argon closes 1 µm below the top, leaving ~0.03 kg (EFCu 231.54 → 176.05 kg).
5. `geom/RTlidfix/` validates the fix (requested by Luigi Pertoldi): `before/` and `after/` hold the GDML built at
   the PR's base and head, `confine.mac` samples 1000 remage vertices uniformly in `reentrance_tube_copper`, and
   `pr142.C` plots the tube's top from both. 240 of 1000 vertices fall in the lid before, none after.

![PR 142](geom/RTlidfix/pr142.png)

```bash
cd geom/RTlidfix
remage -q -g before/l1000.gdml -o before/vtx.root -w confine.mac
remage -q -g after/l1000.gdml -o after/vtx.root -w confine.mac
root -l -b -q pr142.C
```

### 2. Simulation

- **What:** `sim/sim.mac` (remage) decays one nuclide uniformly through one section of the tube wall and records the
  energy each decay leaves in each detector and in the UGLAr.
- **Nuclides:** **Tl-208** (²³²Th chain) and **Bi-214** (²³⁸U chain), the two decays with gammas above `Qββ` that
  are "of primary concern" (`LEGEND1000pCDR, LEGEND.pdf`, p. 29).
- **`-s` aliases:** `GDML`; `Z` `A` (Tl-208: 81 208, Bi-214: 83 214), with `ZMIN` `AMIN` the same so only that
  nuclide decays; `NEV`; `SEED`, which starts remage's random numbers and names the run's file (the same seed repeats
  the same decays, so every run has its own); `VOLS`, a regex over the source volumes `reentrance_tube_copper`
  (EFCu), `reentrance_tube_layer_copper_ofhc` (Cu) or `reentrance_tube_layer_steel_316L` (SS).
- **A run:** 10 M decays of one nuclide in one section, a 32-thread NERSC task of 16-75 min.
- **Output:** `output/output_lh5/<nuclide>_<SEED>.lh5`, kept because every rerun starts here: `stp/` (the steps in
  the germanium detectors and in `liquid_argon_underground`), `vtx` (one row per decay, with its position),
  `detector_origins`, `number_of_simulated_events`.
- **Lid:** no decay is in it. `sim/sim.mac` excludes it; the first six EFCu runs, made before that, had its decays
  removed (7.6 M of each 10 M kept).
- **Why one section at a time:** a run over the whole tube would put nearly 9 in 10 of its decays in the Cu and SS,
  which the detectors barely see. So each section gets the decays it needs, and the reweighting (step 5) combines
  them into any design:

| section | runs | decays (Tl-208 / Bi-214) | real years they stand for | why that many                                      |
| ------- | ---- | ------------------------ | ------------------------- | -------------------------------------------------- |
| EFCu    | 25   | 122 M / 112 M            | 166 000 / 107 000         | 122 of the BI's 128 counts: they set its ±9%       |
| Cu      | 366  | 460 M / 3.2 G            | 55 000 / 116 000          | they set how far the SS can come down (step 5)     |
| SS      | 36   | 50 M / 310 M             | 7 / 6                     | 2.7-4.8 m above the array: no ROI count at any cut |
| total   | 427  | 4.3 × 10⁹ (76 GB)        |                           |                                                    |

- **Real years:** how long the section, as built, takes to make that many decays (its activity × mass).
- **Bi-214 gets ~7 × Tl-208's decays in the Cu and SS:** the SS holds ~7 × more Bi-214 than Tl-208 (2500 vs
  360 µBq/kg), and a Cu count turned into SS by the reweighting must weigh the same in both chains.
- **Cost:** 142 task-hours on NERSC's Perlmutter, 20 node-hours of the m2676 allocation (a Cu run 16 min, an EFCu
  run about an hour); ~220 hours on an M1 Pro MacBook (8 threads). Another 55 node-hours went into runs that did not
  change the answer, and were deleted.

**How many decays each section needs** was measured on the first 408 runs, not guessed:

1. p = simulated counts after all cuts ÷ simulated decays, per nuclide and section. The geometry and the cuts set
   it, not the number of decays, so one or two runs measure it.
2. A nuclide-section that makes a share f of the BI with n counts adds f²/n to the square of the BI's relative
   uncertainty. ±10% needs Σ f²/n ≤ 0.01.
3. Then (A/E cut at −1.06): EFCu 58% with 15 counts (Tl-208) and 33% with 6 (Bi-214), Cu 6.5% with 3 and 2% with 2,
   SS none: Σ f²/n = 0.042, ±21%.
4. Needed: ~80 and ~28 EFCu counts (Σ f²/n = 0.0098), ~5 × the EFCu's decays; the Cu and SS none.
5. Done: 19 more EFCu runs without the lid (10 Tl-208, 9 Bi-214, 2.7 node-hours) brought the EFCu to 122 counts
   after all cuts (77 Tl-208, 45 Bi-214, with the production's A/E cut of −1.5) and the BI to ±9%.
6. Not done: more Cu decays. One ROI count after M1 moved from Cu into SS adds 12% of the BI; bringing that to ~3%
   would take ~4 × the Cu decays (~40 node-hours). Not worth it, so the SS limit stays set by a single count.

To rerun them on Perlmutter (`nersc`: an ssh host alias with an `sshproxy` key; `<SEEDS>`: EFCu 7001-7013 /
8001-8012, Cu 5001-5046 / 6001-6320, SS 3001-3005 / 4001-4031 for Tl-208 / Bi-214; `<NUCLIDE>` `<Z>` `<A>`: `tl208 81 208` or `bi214 83 214`; `<VOLS>` the section's volume):

```bash
ssh nersc 'mkdir -p $SCRATCH/LEGEND1000-Simulation/sim $SCRATCH/LEGEND1000-Simulation/geom/l1000geom $SCRATCH/LEGEND1000-Simulation/output/output_lh5'
rsync -a geom/l1000geom/l1000.gdml nersc:'$SCRATCH/LEGEND1000-Simulation/geom/l1000geom/'
rsync -a sim/sim.mac nersc:'$SCRATCH/LEGEND1000-Simulation/sim/'
ssh nersc 'shifterimg pull docker:legendexp/remage:v1.1.0'
ssh nersc 'cd $SCRATCH/LEGEND1000-Simulation && sbatch --account=m2676 --qos=shared --constraint=cpu --cpus-per-task=32 --time=01:30:00 --output=output/%x_%a.log --image=docker:legendexp/remage:v1.1.0 --array=<SEEDS> --job-name=rt --wrap "shifter --module=none remage -q --ignore-warnings -t 32 -w -m -s GDML=geom/l1000geom/l1000.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/output_lh5/<NUCLIDE>_\$SLURM_ARRAY_TASK_ID.lh5 -s Z=<Z> -s A=<A> -s ZMIN=<Z> -s AMIN=<A> -s VOLS=<VOLS> -- sim/sim.mac"'
rsync -a --include='*_[0-9][0-9][0-9][0-9].lh5' --exclude='*' nersc:'$SCRATCH/LEGEND1000-Simulation/output/output_lh5/' output/output_lh5/
ssh nersc 'rm -r $SCRATCH/LEGEND1000-Simulation'
```

### 3. Detector response

- **What:** `sim/hit.py` ([reboost](https://reboost.readthedocs.io) 1.4.0) turns each run's steps into what is
  measured: in the germanium, the energy through the dead layer, the resolution and the A/E classifier; in the
  UGLAr, the photoelectrons (PE) the SiPMs see.
- **Parameters:** [legend1000-metadata](https://github.com/legend-exp/legend1000-metadata), `simprod/config`,
  experiment `l1000dsg01`: the settings of LEGEND-1000's simulation production
  ([legend-simflow](https://github.com/legend-exp/legend-simflow)).
- **Maps:** Moritz Neuberger's UGLAr optical map v0.4.0 and drift-time maps of `V00000A` (3500 V, ⟨100⟩ and ⟨110⟩
  axes), on NERSC (paths in the commands below).
- **Detectors:** every one is the dummy `V99999Z` with the drift times of `V00000A` (the same geometry); real
  detectors would each bring their own dead layer, resolution and A/E.

|            | value                                                                                 | from                                                         |
| ---------- | ------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| dead layer | full-charge-collection depth (FCCD) 0.7 mm, outer half fully dead, linear in between  | `V99999Z`; `tier/hit` `dead_layer_fraction` 0.5              |
| resolution | FWHM = √(0.5 + 0.001·E) keV: 1.59 keV at `Qββ`                                        | `pars/geds/eresmod`                                          |
| A/E        | A: max of one template current pulse per step, at its drift time; σ_A/E = 0.01        | `pars/geds/currmod`, `aoeresmod`; `output/dtmap_V00000A.lh5` |
| PE         | LAr scintillation × the summed optical map; σ 0.3 PE, 16 ns merging, ≤ 100 PE per hit | `/all` of the optical map; `tier/opt`                        |

- **A/E** is the amplitude of the current pulse over the energy: ≈ 1 for the single-site signature of 0νββ, << 1
  for a gamma that interacts several times (`LEGEND1000pCDR, LEGEND.pdf`, p. 24). One template pulse is used, at
  each step's drift time, not a pulse-shape library.
- **Optical map:** for every 5 mm cell of the UGLAr, the probability that a scintillation photon made there reaches
  a SiPM; UGLAr scintillation only (`Optical_map_L1000, Sanchez.pdf`, p. 2). Its `/all` group sums every SiPM, as
  the LAr veto does (the production's second condition, 4 SiPMs with light, and random coincidences are left out).
- **The map's limits:** made for legend-pygeom-l1000 0.4.0, it stops 68 mm short of this UGLAr's outer wall and
  gives no light above z = 1.375 m (3.61 m deep). Until a 0.6.0 map arrives, the BI is a 0.9-2.8 × 10⁻⁸ range
  (2.8 × 10⁻⁸ with no light from the argon beyond the map's edge).
- **Output:** `output/output_root/<nuclide>_<SEED>_hit.root` (~7 s a file, 20-40 s for an EFCu run; ~50 GB in all;
  only `ana/ana.C` reads them): TTrees `geds` (`evtid`, `det`, `edep`, `energy` [keV], `aoe_class`), `lar` (`evtid`,
  `edep`, `pe`, `pe_map`), `vtx` (`evtid`, `xloc`, `yloc`, `zloc` [m]).
- **Options:** `--optmap`, `--dtmap`, `--metadata`, `--gdml`, `--outdir` (defaults: the two maps fetched below, the
  metadata above, `geom/l1000geom/l1000.gdml`, `output/output_root`); `--seed` (default 1; not remage's `SEED`: with
  each file's name it seeds the smearing, photons and photoelectrons, so a file always gives the same hits);
  `--force` (rebuilds hit files that are already newer than their inputs; without it they are skipped).
- **A new optical map** needs only `--optmap` and the analysis (step 4) again, no new simulation.

The two maps, then the response:

```bash
rsync -a --partial nersc:/global/cfs/cdirs/m2676/users/neuberger/L1000_optical_muon_sims/omaps/v0.4.0/ular/merged/merged_optmap_20260225_063750.lh5 output/
rsync -a nersc:/global/cfs/cdirs/m2676/users/neuberger/L1000_optical_muon_sims/hpge_related/dtmaps/gen/V00000A.lh5 output/dtmap_V00000A.lh5
~/venvs/v/bin/python sim/hit.py output/output_lh5/*.lh5
```

### 4. Cuts and BI

`ana/ana.C` is the analysis and its figures:

- **Arguments,** both optional: the hit files, as `pattern=nuclide` pairs (wildcards merge files; default: the Tl-208
  and Bi-214 sets in `output/output_root`), and the tube table (default `output/tube.txt`, from step 1).
- **Settings:** the cuts, activities and the 10% line are in `cfg` (`m1_keV`, `lar_pe`, `psd_low`, `chain`, `simTol`).

**Cuts**, in the order of LEGEND-1000's background model (`CDRBackgroundModel, Green.pdf`, p. 6), with the
production's values (legend1000-metadata, `tier/evt` and `pars/geds/psdcuts`):

1. **M1:** exactly one detector above 25 keV.
2. **LAr veto:** rejects ≥ 4 PE.
3. **A/E** > −1.5, single-site-like: one value for all detectors.

**Is −1.5 the right A/E cut?**

1. In data, LEGEND tunes the cut to keep 90% of the Tl-208 double-escape peak (DEP, 1592.5 keV, single-site)
   (`LEGEND1000pCDR, LEGEND.pdf`, p. 71).
2. Here −1.5 keeps 87 ± 3% of the 2992 simulated decays in the DEP after M1 (net of the continuum).
3. Keeping 90% would take −6.6, because this one-template A/E leaves the DEP's A/E with a long flat tail.
4. At −6.6 the other peaks pass 3 × what was measured for `V00000A` at HADES, the underground lab where LEGEND
   characterises its detectors (the survival fractions in Moritz Neuberger's `definition/diodes/V00000A.yaml`).
5. So the production's cut is kept. What it keeps of each Tl-208 peak and of Tl-208's ROI counts, after M1:

| kept by A/E        | DEP 1592.5 keV | SEP 2103.5 keV (single escape) | FEP 2614.5 keV (full energy) | ROI  |
| ------------------ | -------------- | ------------------------------ | ---------------------------- | ---- |
| this simulation    | 0.87           | 0.03                           | 0.05                         | 0.29 |
| `V00000A` at HADES | 0.90           | 0.04                           | 0.06                         | 0.27 |

**BI:**

- **ROI:** one window, `Qββ` ± 50 keV (1989-2089 keV, **100 keV**), stopping short of Tl-208's single-escape peak
  at 2103.5 keV. Ralph Massarczyk's tube study uses ± 55 keV (`ralph_GeneralCall_10_12_2026.pdf`, 13 Oct 2026, p. 3).
- **BI** = counts in the ROI per keV, per kg of germanium (1023.4 kg), per year: counts/(keV·kg·yr).
- **Goal:** LEGEND-1000's is below 10⁻⁵ counts/(keV·kg·yr) (`LEGEND1000pCDR, LEGEND.pdf`, p. 8): 1 count a year in
  10³ kg and 10² keV.
- **Signal window:** `Qββ` ± 2σ (2037.6-2040.4 keV); the background expected there is the BI × 2.7 keV × the exposure.

**From decays to counts.** Each simulated decay stands for w real decays a year, its weight (the same weight as in
`ralph_GeneralCall_10_12_2026.pdf`, p. 3):

- w = real decays / simulated decays = mass (kg) × activity (µBq/kg) × 0.36 (Tl-208) or 1 (Bi-214) × 1 yr /
  simulated decays, knowing Bq = decays/s. Tl-208 is 35.94% of the ²³²Th chain
  (`Nuclear Data Sheets for A = 208, Martin.pdf`, p. 1); Bi-214 all of the ²³⁸U chain;
- Σw over the simulated decays that land in the ROI = the counts a year;
- BI = counts a year / (1023.4 kg × 100 keV).

Bulk activities only, chains in equilibrium; the mass of each section is the GDML's. The activities (`cfg::chain`):

| µBq/kg | ²³²Th chain | ²³⁸U chain | source                                                                                                                                                    |
| ------ | ----------- | ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| SS     | 1000        | 2500       | Bernhard Schwingenheuer, quoted in `RT MaterialMix, ralph.pdf` (23 Jun 2026), p. 5                                                                        |
| Cu     | 1.1         | 1.3        | MAJORANA's assay, as quoted there (p. 15); the paper itself gives 1.10 ± 0.14 and 1.37 ± 0.18 (`An assay-based background projection, Majorana Demonstrator.pdf`, p. 6) |
| EFCu   | 0.37        | 0.19       | CDR Table XVII (`CDRBackgroundModel, Green.pdf`, p. 7)                                                                                                    |

The tube as built, from mass to BI (each section sampled at its own density; decays counted where they happened,
so ~1 M of the EFCu runs' decays, just above the Cu seam, count as Cu):

| chain  | section | mass     | real decays a year | simulated decays | w          | simulated counts: no cut / M1 / + LAr / + A/E | counts a year, all cuts | BI, all cuts       |
| ------ | ------- | -------- | ------------------ | ---------------- | ---------- | --------------------------------------------- | ----------------------- | ------------------ |
| Tl-208 | EFCu    | 174.6 kg | 733                | 122 M            | 6.0 × 10⁻⁶ | 63755 / 52911 / 437 / 77                      | 4.6 × 10⁻⁴              | 4.52 ± 0.52 × 10⁻⁹ |
| Tl-208 | Cu      | 671.7 kg | 8 380              | 461 M            | 1.8 × 10⁻⁵ | 551 / 497 / 34 / 3                            | 0.5 × 10⁻⁴              | 0.53 ± 0.31 × 10⁻⁹ |
| Tl-208 | SS      | 619.4 kg | 7.0 × 10⁶          | 50.0 M           | 0.14       | 0 / 0 / 0 / 0                                 | 0                       | 0                  |
| Bi-214 | EFCu    | 174.6 kg | 1 047              | 112 M            | 9.4 × 10⁻⁶ | 4720 / 4383 / 588 / 45                        | 4.2 × 10⁻⁴              | 4.11 ± 0.61 × 10⁻⁹ |
| Bi-214 | Cu      | 671.7 kg | 27 560             | 3.20 G           | 8.6 × 10⁻⁶ | 170 / 160 / 42 / 3                            | 0.3 × 10⁻⁴              | 0.25 ± 0.15 × 10⁻⁹ |
| Bi-214 | SS      | 619.4 kg | 4.9 × 10⁷          | 310 M            | 0.16       | 0 / 0 / 0 / 0                                 | 0                       | 0                  |
| all    |         |          |                    | 4.26 G           |            | 128 after all cuts                            | 9.6 × 10⁻⁴              | 9.4 ± 0.9 × 10⁻⁹   |

```bash
root -l -b -q ana/ana.C
```

### 5. Reweighting

**A design** moves the two seams:

- the SS ends at depth L1 (2.06 m as built), the EFCu starts at depth L2 (4.07 m), with Cu between;
- the wall keeps its shape and thickness as built; only the material of each part changes (EFCu turned into Cu at
  4.07-4.40 m is the EFCu's 1.5 mm wall, made of Cu);
- `ana/ana.C` scans them: for every L2, L1 goes down 1 cm at a time while the BI stays within 10% of the as-built's;
- the answer is the least EFCu, then, for it, the most SS.

**Reweighting, not resimulating:**

1. What a decay does in the detectors depends on where it happens, and very little on what the wall there is made
   of: 5 mm of wall absorbs 14% of 2.6 MeV gammas in steel and 16% in Cu (`ralph_GeneralCall_10_12_2026.pdf`, p. 5).
2. So each simulated decay can stand for any material a design puts at its depth; only its weight w changes.
3. The new w = that material's real decays a year per m³ (density × activity) ÷ the simulated decays per m³ of the
   run the decay came from.
4. A design's counts a year are Σw with the new weights, so any L1 and L2 can be tried without a new simulation.
5. The decays at a depth come from the section built there (SS runs 0-2.06 m, Cu runs 2.06-4.07 m, EFCu runs
   4.07-6.25 m), so the two savings reweight different runs:

| design change                     | decays reweighted                                  | w × (Tl-208 / Bi-214) | simulated counts there           |
| --------------------------------- | -------------------------------------------------- | --------------------- | -------------------------------- |
| EFCu starts lower (L2 deeper)     | the EFCu runs', from 4.07 m to L2: EFCu becomes Cu | 3.0 / 6.8             | thousands: the change is precise |
| SS comes further down (L1 deeper) | the Cu runs', from 2.06 m to L1: Cu becomes SS     | 810 / 1700            | a handful: each one matters      |

**Why the Cu runs decide the SS:**

1. SS added in a design sits where the Cu is now, 2-3 m above the array; only the Cu runs have decays there.
2. Those decays rarely reach the detectors: one in ~5 × 10⁷ Tl-208 decays gives an ROI count after M1, and none of
   1.7 × 10⁹ Bi-214 decays does (the table below).
3. Turned into SS, each such count weighs ~1000 times more: one adds 12% of the BI, more than the 10% line.
4. So the SS can come down to 2.69 m, just above the first such count (2.70 m): the step in the design map.

**Cut survival:**

1. Away from the detectors, too few counts survive all cuts to measure a design's change on them.
2. So the change is measured on the ROI counts that pass M1 (100-700 times more) × the **cut survival**: the
   fraction of them that also pass the LAr veto and A/E, per 0.62 m slab (the table):
   - near the detectors (below 4.37 m), per chain: 1.4 × 10⁻³ (Tl-208), 1.0 × 10⁻² (Bi-214);
   - higher up, both chains together: 10 of 1258 = 8.0 × 10⁻³ (~4 × Tl-208's value near the detectors, partly
     because the optical map gives no light above 3.61 m deep).

| depth [m]                        | Tl-208 decays | ROI after M1 | after all cuts | Bi-214 decays | ROI after M1 | after all cuts |
| -------------------------------- | ------------- | ------------ | -------------- | ------------- | ------------ | -------------- |
| 0.00-1.87                        | 45.4 M        | 0            | 0              | 282 M         | 0            | 0              |
| 1.87-2.50                        | 105 M         | 0            | 0              | 729 M         | 0            | 0              |
| 2.50-3.12                        | 144 M         | 5            | 0              | 1.00 G        | 0            | 0              |
| 3.12-3.75                        | 144 M         | 97           | 0              | 1.00 G        | 29           | 2              |
| 3.75-4.37                        | 90.3 M        | 959          | 6              | 511 M         | 168          | 2              |
| 4.37-5.00                        | 36.7 M        | 11645        | 16             | 33.8 M        | 958          | 7              |
| 5.00-5.62                        | 36.7 M        | 26443        | 37             | 33.7 M        | 2217         | 23             |
| 5.62-6.25, nearest the detectors | 29.9 M        | 14259        | 21             | 27.5 M        | 1171         | 14             |

## Results

**BI as built** [counts/(keV·kg·yr)]

| no cut    | M1                         | M1 + LAr                         | all cuts (M1 + LAr + A/E)       |
| --------- | -------------------------- | -------------------------------- | ------------------------------- |
| 4.29×10⁻⁶ | 3.61×10⁻⁶ (M1 removes 16%) | 0.89×10⁻⁷ (LAr veto removes 98%) | **0.94×10⁻⁸** (A/E removes 89%) |

The tube's spectrum as built, through each cut: the ROI shaded green with `Qββ` as a line, the BI after all cuts as
the black line across it, the 10⁻⁵ goal dashed.

![the spectrum](output/slide_spectrum.png)

**Designs**

For each L2, the deepest L1 that keeps the BI within 10%; deeper than 4.40 m no L1 passes:

| EFCu from (L2)         | SS to (L1) | SS : Cu : EFCu          | BI change |
| ---------------------- | ---------- | ----------------------- | --------- |
| 4.07 m, as built       | 2.69 m     | 809 : 460 : 175 kg      | none      |
| 4.20 m                 | same       | same : 472 : 163 kg     | +2%       |
| 4.30 m                 | same       | same : 480 : 154 kg     | +4%       |
| **4.40 m** (−17% EFCu) | same       | **same : 489 : 146 kg** | **+8%**   |

![the designs](output/slide_designs.png)

The rt viewer of `/Applications/LEGEND.app` shows the same scan in 3D; rebuild it after each `ana/ana.C` run:

```bash
~/venvs/v/bin/python app/make_app.py
```

## LEGEND 1700?

Matthew Busch's talk, _String and lock updates, including a 12 detector string concept for beyond 1T deployment_
(29 Sep 2026; `Hall C design, Beusch.pdf`), asks, at Bernhard's request: could the same hardware later hold more
than 1 tonne of germanium?

- **The idea** (p. 2): 12 detectors per string instead of 8, at the same central depth in the argon: 42 × 8 × 3 kg =
  1008 kg becomes 48 × 12 × 3 kg ≈ 1.7 tonnes.
- **What changes** (p. 3, 8, 14): strings get 0.28 m longer at each end, the tube gets 0.2 m more OFHC copper and a
  deeper bottom head (467 → 575 mm), and the cleanroom floor, cryostat and water tank come down 0.2 m.
- **Status** (p. 14): viable, but still waiting on the crane and cleanroom details.
- **The catch** (p. 14): he recommends the Hall C change only if it works with the background simulations.

A possible next step for this study is that check: build the 12-detector geometry in legend-pygeom-l1000, then rerun
the same simulation and analysis on it.
