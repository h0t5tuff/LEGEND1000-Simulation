# LEGEND-1000 re-entrance tube

The re-entrance tube (RT) holds the underground argon (UGLAr) and the detector array. It is a 6 m × 2 m cylinder of
steel (SS) at the top, oxygen-free high-conductivity (OFHC) copper (Cu) in the middle and electroformed copper (EFCu),
3-7 × less radioactive than the Cu, around the detectors.

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

What EFCu's own radioactivity looks like in the detectors. The 1008 EFCu detector holders
(`hpge_string_support_weldment_copper`, 1.92 kg in all), clamping each detector's PEN plate to the support rods, decay with
EFCu's Th-232 and U-238 chains: 370 k and 190 k chains, the ratio of EFCu's activities (0.37 : 0.19 µBq/kg), so the
two spectra add as counted. Two lines threaten `Qββ` = 2039 keV, the energy of neutrinoless double-beta decay
(0νββ): **Tl-208 at 2614 keV** and **Bi-214 at 2448 keV**, whose gammas can leave that much in a detector.

![the holders](output/holders_geometry.png)

![the holders' spectrum](output/holders_spectrum.png)

The energy each of the 336 detectors gets, without the detector response (§3), summed over 10 µs: 434 263 such
detector events from the Th-232 chains, 176 339 from the U-238. `sim/holders.mac` writes every germanium step (`stp/germanium`: `evtid`, `det_uid`, `edep_in_keV`,
`time_in_ns`, ...) and one row per chain (`stp/vtx`) to `output/holders_th232.root` and `holders_u238.root`;
`ana/holders.C` draws them.

```bash
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --holders output/holders_geometry.png
remage -q --ignore-warnings -t 8 -w -o output/holders_th232.root -s GDML=geom/l1000geom/l1000.gdml -s Z=90 -s A=232 -s NEV=370000 -s SEED=21 -- sim/holders.mac
remage -q --ignore-warnings -t 8 -w -o output/holders_u238.root -s GDML=geom/l1000geom/l1000.gdml -s Z=92 -s A=238 -s NEV=190000 -s SEED=22 -- sim/holders.mac
root -l -b -q ana/holders.C
```

## Method

Five steps. From the NERSC decays already in `output/output_lh5`, the commands of §3 and §4 redo the response and
the analysis, §5 included (it runs inside `ana/background.C`).

### 1. The geometry

`geom/l1000geom/l1000.gdml` (legend-pygeom-l1000 0.6.0), the full LEGEND-1000 geometry. The tube is 6.246 m tall and
0.999 m in radius, its walls EFCu 1.5 mm (2.1 mm at the Cu seam), Cu and SS 6 mm; the 336 detectors sit 4.8-5.9 m
below its top, inside the EFCu section. `ana/background.C` reads all of this from the GDML, with the densities and
the germanium mass, and checks it before using it: a clean outline, the bottom head and the Cu and SS layers flush
with the barrel.

![the tube](output/tube.png)

The tube's outline as `ana/background.C` reads it: its three sections, the decay sources, and the two seams a design
moves (z = 920 and 2925 mm in the GDML).

![the tube, the decay source](output/tube_geometry.png)

![the UGLAr, the veto](output/uglar_geometry.png)

The underground argon around the detectors (18.83 m³, 26.17 t) is the veto, its light read by 12096 fibres and 252
silicon photomultipliers (SiPMs).

![the HPGe detectors](output/hpge_geometry.png)

Each detector is an inverted-coaxial point-contact (ICPC) germanium detector, 90.0 mm tall and 88.8 mm across, its
borehole 53 mm deep from the top: 3.046 kg, 1023.4 kg in all. `geom/l1000geom/l1000-autopeel.py` draws each part
alone, offscreen, and prints its numbers:

```bash
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --tube output/tube_geometry.png
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --uglar output/uglar_geometry.png
~/venvs/v/bin/python geom/l1000geom/l1000-autopeel.py --hpge output/hpge_geometry.png
```

### 2. Simulation

With `sim/run.mac` (remage), one nuclide decays uniformly through one section of the tube wall: **Tl-208** (²³²Th
chain) or **Bi-214** (²³⁸U chain), the only chain members whose gammas reach `Qββ`. It records the energy each decay
leaves in each detector and in the UGLAr. Its `-s` aliases: `GDML`; `Z` `A` (Tl-208: 81 208, Bi-214: 83 214); `NEV`;
`SEED`, which starts remage's random numbers (the same seed repeats the same decays, so every run has its own) and
names the run's file; `VOLS`, a regex over the source volumes `reentrance_tube_copper` (EFCu),
`reentrance_tube_layer_copper_ofhc` (Cu) or `reentrance_tube_layer_steel_316L` (SS). Each run is 10 M decays (7.6 M
in the EFCu runs) of one nuclide in one section, a 32-thread NERSC task of 16-60 min, and writes
`output/output_lh5/<nuclide>_<SEED>.lh5`: `stp/` (the steps in the germanium detectors and in
`liquid_argon_underground`), `vtx` (one row per decay, with its position), `detector_origins`,
`number_of_simulated_events`. They are kept: every rerun starts here.

The sections are simulated separately because their decays reach the detectors at very different rates, so each gets
the decays it needs (a run over the whole tube would put nearly 9 in 10 of its decays in the Cu and SS, which the detectors
barely see); §5 combines them into any design:

| section | runs | decays (Tl-208 / Bi-214) | real years they stand for | why that many                                                                                           |
| ------- | ---- | ------------------------ | ------------------------- | ------------------------------------------------------------------------------------------------------- |
| EFCu    | 6    | 23 M / 23 M              | 31 000 / 22 000           | next to the array, they give 94 of the BI's 113 events after all cuts: the BI's ±10%                    |
| Cu      | 366  | 460 M / 3.2 G            | 55 000 / 116 000          | they decide how far the SS can come down (§5): enough that one event, turned into SS, adds 3% of the BI |
| SS      | 36   | 50 M / 310 M             | 7 / 6                     | 2.7-4.8 m above the array: none of their decays puts an event in the ROI (§4), at any cut level         |
| total   | 408  | 4.1 × 10⁹                |                           | 66 GB                                                                                                   |

"Real years": how long the section, as built, takes to make that many decays (its activity × mass, §4). In the Cu and
SS, Bi-214 needs ~7 × the decays of Tl-208: the SS holds ~7 × more Bi-214 than Tl-208 (2500 vs 360 µBq/kg), and an
event turned into SS (§5) must weigh the same in both chains. On NERSC's Perlmutter the runs took 122 task-hours, 17
node-hours of the m2676 allocation (16 min for a typical Cu run); back to back on an M1 Pro MacBook (8 threads,
22.5 min for that run) ~175 hours, about a week. Another 55 node-hours went into runs that turned out not to change
the answer, and were deleted.

To rerun them (`<SEEDS>`: EFCu 7001-7003 / 8001-8003, Cu 5001-5046 / 6001-6320, SS 3001-3005 / 4001-4031 for
Tl-208 / Bi-214; `<NUCLIDE>` `<Z>` `<A>`: `tl208 81 208` or `bi214 83 214`; `<VOLS>` the section's volume):

```bash
ssh nersc 'mkdir -p $SCRATCH/LEGEND1000-Simulation/sim $SCRATCH/LEGEND1000-Simulation/geom/l1000geom $SCRATCH/LEGEND1000-Simulation/output/output_lh5'
rsync -a geom/l1000geom/l1000.gdml nersc:'$SCRATCH/LEGEND1000-Simulation/geom/l1000geom/'
rsync -a sim/run.mac nersc:'$SCRATCH/LEGEND1000-Simulation/sim/'
ssh nersc 'shifterimg pull docker:legendexp/remage:v1.1.0'
ssh nersc 'cd $SCRATCH/LEGEND1000-Simulation && sbatch --account=m2676 --qos=shared --constraint=cpu --cpus-per-task=32 --time=01:30:00 --output=output/%x_%a.log --image=docker:legendexp/remage:v1.1.0 --array=<SEEDS> --job-name=rt --wrap "shifter --module=none remage -q --ignore-warnings -t 32 -w -m -s GDML=geom/l1000geom/l1000.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/output_lh5/<NUCLIDE>_\$SLURM_ARRAY_TASK_ID.lh5 -s Z=<Z> -s A=<A> -s VOLS=<VOLS> -- sim/run.mac"'
rsync -a --include='*_[0-9][0-9][0-9][0-9].lh5' --exclude='*' nersc:'$SCRATCH/LEGEND1000-Simulation/output/output_lh5/' output/output_lh5/
ssh nersc 'rm -r $SCRATCH/LEGEND1000-Simulation'
```

### 3. Detector response

`sim/hit.py` (reboost 1.4) turns each run's steps into what the detectors measure: in the germanium, the energy
through the dead layer, the resolution and the A/E classifier; in the UGLAr, the photoelectrons (PE) the SiPMs see.
Its inputs: the response parameters of [legend1000-metadata](https://github.com/legend-exp/legend1000-metadata)
(`~/Documents/legend-exp/legend1000-metadata`, `simprod/config`, experiment `l1000dsg01`: the settings of
LEGEND-1000's simulation production, which [legend-simflow](https://github.com/legend-exp/legend-simflow) runs), the UGLAr optical map v0.4.0 and the drift-time maps
of `V00000A` (3500 V, ⟨100⟩ and ⟨110⟩ axes), applied with [reboost](https://reboost.readthedocs.io). Every detector is
the dummy `V99999Z`, with the drift times of `V00000A` (the same geometry): a stand-in, since real detectors would each
bring their own dead layer, resolution and A/E.

|            | value                                                                                 | from                                                         |
| ---------- | ------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| dead layer | full-charge-collection depth (FCCD) 0.7 mm, outer half fully dead, linear in between  | `V99999Z`; `tier/hit` `dead_layer_fraction` 0.5              |
| resolution | FWHM = √(0.5 + 0.001·E) keV: 1.59 keV at `Qββ`                                        | `pars/geds/eresmod`                                          |
| A/E        | A: max of one template current pulse per step, at its drift time; σ_A/E = 0.01        | `pars/geds/currmod`, `aoeresmod`; `output/dtmap_V00000A.lh5` |
| PE         | LAr scintillation × the summed optical map; σ 0.3 PE, 16 ns merging, ≤ 100 PE per hit | `/all` of the optical map; `tier/opt`                        |

- **A/E** is the amplitude of the current pulse over the energy. A gamma that scatters several times in the crystal
  spreads its charge out in time, giving a broad, low pulse: A/E << 1. 0νββ puts its energy (two electrons) in one
  spot: A/E ≈ 1. The drift-time maps give how long each step's charge takes to reach the contact, which shapes the
  pulse; one template pulse is used, not a pulse-shape library.
- **Optical map** (Moritz Neuberger's, on NERSC): for every 5 mm cell of the UGLAr, the probability that a
  scintillation photon made there reaches a SiPM; its `/all` group sums every SiPM, as the LAr veto does (the
  production's second condition, 4 SiPMs with light, and random coincidences are left out). This v0.4.0 map stops
  68 mm short of this UGLAr's outer wall (the argon there takes the light at the map's edge), and argon above
  z = 1.375 m gives no light, so until the 0.6.0 map
  arrives the BI is a 1.0-3.2 × 10⁻⁸ range (3.2 × 10⁻⁸ with the LAr veto read off the map alone, no light from the
  argon beyond its edge).

`sim/hit.py` writes `output/output_root/<nuclide>_<SEED>_hit.root` (~7 s a file, ~50 GB in all; only `ana/background.C` reads
them): TTrees `geds` (`evtid`, `det`, `edep`, `energy` [keV], `aoe_class`), `lar` (`evtid`, `edep`, `pe`, `pe_map`),
`vtx` (`evtid`, `xloc`, `yloc`, `zloc` [m]). Options: `--optmap`, `--dtmap`, `--metadata`, `--gdml`, `--outdir`
(defaults: the two maps fetched below, the metadata above, §1's GDML, `output/output_root`); `--seed` (default 1; not remage's `SEED`: with each
file's name it seeds the smearing, photons and photoelectrons, so a file always gives the same hits); `--force`
(rebuilds a hit file that is newer than its inputs, otherwise kept). A new optical map needs only `--optmap` and §4
again, no new simulation. The two maps, then the response:

```bash
rsync -a --partial nersc:/global/cfs/cdirs/m2676/users/neuberger/L1000_optical_muon_sims/omaps/v0.4.0/ular/merged/merged_optmap_20260225_063750.lh5 output/
rsync -a nersc:/global/cfs/cdirs/m2676/users/neuberger/L1000_optical_muon_sims/hpge_related/dtmaps/gen/V00000A.lh5 output/dtmap_V00000A.lh5
~/venvs/v/bin/python sim/hit.py output/output_lh5/*.lh5
```

### 4. Cuts and BI

`ana/background.C` is the analysis and its figures. Its arguments, both optional: the hit files, as `pattern=nuclide`
pairs (wildcards merge files; default: the Tl-208 and Bi-214 sets in `output/output_root`), and the GDML (default
`geom/l1000geom/l1000.gdml`). The cuts, activities and the 10% line are in `cfg` (`m1_keV`, `lar_pe`, `depKeep`,
`chain`, `simTol`).

**Cuts**, in LEGEND's order: **M1** (exactly one detector above 25 keV) and the **LAr veto** (rejects ≥ 4 PE), both
from legend1000-metadata's `tier/evt`; then **A/E** > −1.06, single-site-like, placed by LEGEND's convention to keep
90% of the Tl-208 double-escape peak (DEP, 1592.5 keV, single-site; the production's −1.5 keeps 92%). It is one value
for all detectors, tuned here on the simulation's 488 DEP events after M1 (net of the continuum; the 90% is ±1%
statistics). What it keeps of each Tl-208 peak and of the ROI, next to `V00000A`'s measurement at HADES (the
underground lab where LEGEND characterises its detectors):

| kept by A/E        | DEP 1592.5 keV | SEP 2103.5 keV (single escape) | FEP 2614.5 keV (full energy) | ROI  |
| ------------------ | -------------- | ------------------------------ | ---------------------------- | ---- |
| this simulation    | 0.90           | 0.02                           | 0.04                         | 0.25 |
| `V00000A` at HADES | 0.90           | 0.04                           | 0.06                         | 0.27 |

**BI** is counts in the **ROI**, MAJORANA's background region ([Haufe et al.](https://arxiv.org/abs/2209.10592),
[MAJORANA's final result](https://arxiv.org/abs/2207.07638)): 1950-2350 keV, less ±5 keV around each of 2039.0
(`Qββ`), 2103.5 (Tl-208 single escape), 2118.5 and 2204.1 keV (Bi-214), so **360 keV** of continuum, per keV, per kg
of germanium (1023.4 kg), per year. The wide ROI collects enough events, and the continuum there is flat, so its
average holds at `Qββ`. A signal would sit in `Qββ` ± 2σ (2037.6-2040.4 keV, ~95% of a 0νββ peak at this
resolution): the background expected there is the BI × 2.7 keV × the exposure.

**Normalisation** turns simulated decays into real ones. Per section, from the GDML and the activities: M = ρV,
A = M × a (a: the chain's µBq/kg × the nuclide's share of it), N_real = A × 1 yr. The simulation put N_sim decays
there, so each stands for N_real / N_sim real decays a year: its weight, and an ROI event adds weight / (1023.4 kg ×
360 keV) to the BI. Bulk activities only, chains in equilibrium. The activities (`cfg::chain`) are from Ralph
Massarczyk's MaterialMix talk (23 Jun 2026); Bi-214 is 100% of the ²³⁸U chain, Tl-208 35.94% of the ²³²Th chain
([Nuclear Data Sheets](https://doi.org/10.1016/j.nds.2007.07.001)):

| µBq/kg | ²³²Th chain | ²³⁸U chain | source quoted there                                           |
| ------ | ----------- | ---------- | ------------------------------------------------------------- |
| SS     | 1000        | 2500       | Bernhard                                                      |
| Cu     | 1.1         | 1.3        | [MAJORANA assay](https://doi.org/10.1103/PhysRevC.110.055804) |
| EFCu   | 0.37        | 0.19       | M. Green                                                      |

and the weights they give (each section sampled at its own density, §2):

| section | V [m³]  | M [kg] | simulated decays (Tl-208 / Bi-214) | weight (Tl-208 / Bi-214)  |
| ------- | ------- | ------ | ---------------------------------- | ------------------------- |
| EFCu    | 0.01948 | 174.6  | 22.6 M / 22.6 M                    | 3.24 × 10⁻⁵ / 4.63 × 10⁻⁵ |
| Cu      | 0.07497 | 671.7  | 460 M / 3.20 G                     | 1.82 × 10⁻⁵ / 8.61 × 10⁻⁶ |
| SS      | 0.07742 | 619.4  | 50.0 M / 310 M                     | 0.140 / 0.158             |

```bash
root -l -b -q ana/background.C
```

### 5. Finding the best design: reweighting

A design moves the two seams: the SS ends at depth L1 (2.06 m as built) and the EFCu starts at depth L2 (4.07 m),
with Cu between. The wall keeps its shape and thickness as built; only the material of each part changes (EFCu turned
into Cu at 4.07-4.50 m is the EFCu's 1.5 mm wall, made of Cu). `ana/background.C` scans them: for every L2, L1 goes
down 1 cm at a time while the BI stays within 10% of the as-built's. The answer is the least EFCu, then, for it, the
most SS.

**Reweighting, not resimulating.** What a decay does in the detectors depends on where it happens, and very little on
what the wall there is made of: SS and Cu stop gammas about equally (8.0 vs 9.0 g/cm³). So each simulated decay can
stand for any material a design puts at its depth, with that material's weight (§4): its real decays a year per m³
(density × activity) over the simulated decays per m³ of the run the decay came from. A design's BI is the sum of its
events' weights, so any L1 and L2 can be tried without a new simulation. The decays at a depth come from the section
built there (SS runs 0-2.06 m, Cu runs 2.06-4.07 m, EFCu runs 4.07-6.25 m), so the two savings reweight different
runs:

| design change                     | decays reweighted                                  | weight × (Tl-208 / Bi-214) | events there                     |
| --------------------------------- | -------------------------------------------------- | -------------------------- | -------------------------------- |
| EFCu starts lower (L2 deeper)     | the EFCu runs', from 4.07 m to L2: EFCu becomes Cu | 3.0 / 6.8                  | thousands: the change is precise |
| SS comes further down (L1 deeper) | the Cu runs', from 2.06 m to L1: Cu becomes SS     | 810 / 1700                 | a handful: each one counts       |

**Why the Cu runs decide the SS.** The SS added in a design would sit where the Cu is now, 2-3 m above the array,
and only the Cu runs have decays there. Those decays rarely reach the detectors: one in 10⁷-10⁸ Tl-208 decays and one
in ~10⁹ Bi-214 decays gives an ROI event after M1 (the table below). Turned into SS, each such event counts ~1000 times
more. With 460 M / 3.2 G Cu decays, one of them adds 3.2% of the BI, so the 10% line is about three events: a
measurement, not one chance event. With only a tenth of them, one event would add ~30%; none happened to fall
between 2.06 and 2.66 m, and the SS looked free down to 2.66 m: the wrong answer.

**Cut survival.** Too few events survive all cuts away from the detectors to count the change on them, so it is
counted on the ROI events that pass M1, 100-500 times more, times the **cut survival**: the fraction of them that also
pass the LAr veto and A/E. It is measured where the change is, from the ROI events per 0.62 m slab below the tube's
top (the table):

- near the detectors (below 4.37 m), per chain: 1.9 × 10⁻³ (Tl-208), 8.6 × 10⁻³ (Bi-214);
- higher up, both chains together: 20 of 2559 = 7.8 × 10⁻³. Tl-208's events survive ~4 × more often there, partly
  because the optical map gives no light above z = 1.375 m.

| depth [m]                        | Tl-208 decays | ROI after M1 | after all cuts | Bi-214 decays | ROI after M1 | after all cuts |
| -------------------------------- | ------------- | ------------ | -------------- | ------------- | ------------ | -------------- |
| 0.00-1.87                        | 45.4 M        | 0            | 0              | 282 M         | 0            | 0              |
| 1.87-2.50                        | 105 M         | 1            | 0              | 729 M         | 0            | 0              |
| 2.50-3.12                        | 144 M         | 13           | 0              | 1.00 G        | 2            | 0              |
| 3.12-3.75                        | 144 M         | 322          | 6              | 1.00 G        | 64           | 1              |
| 3.75-4.37                        | 74.5 M        | 1813         | 10             | 496 M         | 344          | 3              |
| 4.37-5.00                        | 6.83 M        | 8796         | 18             | 6.82 M        | 445          | 4              |
| 5.00-5.62                        | 6.82 M        | 20243        | 42             | 6.82 M        | 1077         | 9              |
| 5.62-6.25, nearest the detectors | 5.55 M        | 11024        | 15             | 5.56 M        | 574          | 5              |

## Results

**The BI as built** (the LEGEND-1000 goal is 10⁻⁵), then per chain and section (the SS has no ROI event at any cut
level):

| no cut      | M1                          | M1 + LAr                                     | all cuts (M1 + LAr + A/E)                           |
| ----------- | --------------------------- | -------------------------------------------- | --------------------------------------------------- |
| 4.60 × 10⁻⁶ | 3.92 × 10⁻⁶: M1 removes 15% | 1.11 ± 0.03 × 10⁻⁷: the LAr veto removes 97% | **0.98 ± 0.10 × 10⁻⁸**: A/E removes 91% of the rest |

| chain  | section | ROI events: no cut / M1 / + LAr / + A/E | BI after all cuts  |
| ------ | ------- | --------------------------------------- | ------------------ |
| Tl-208 | EFCu    | 47851 / 40461 / 457 / 76                | 6.68 ± 0.77 × 10⁻⁹ |
| Tl-208 | Cu      | 1937 / 1751 / 176 / 15                  | 0.74 ± 0.19 × 10⁻⁹ |
| Bi-214 | EFCu    | 2270 / 2110 / 470 / 18                  | 2.26 ± 0.53 × 10⁻⁹ |
| Bi-214 | Cu      | 414 / 396 / 131 / 4                     | 0.09 ± 0.05 × 10⁻⁹ |

![the spectrum](output/slide_spectrum.png)

The tube's spectrum as built, through each cut: the ROI shaded light green, `Qββ` ± 2σ dark green, the BI after all
cuts as the black line across the ROI, the 10⁻⁵ goal dashed.

**The designs**: for each depth L2 where the EFCu starts, the deepest L1 where the SS can end:

| EFCu from (L2)    | EFCu       | SS to (L1) | SS : Cu : EFCu         | BI change                      |
| ----------------- | ---------- | ---------- | ---------------------- | ------------------------------ |
| 4.07 m, as built  | 175 kg     | 2.69 m     | 809 : 460 : 175 kg     | +6%                            |
| 4.20 m            | 163 kg     | 2.69 m     | 809 : 472 : 163 kg     | +8%                            |
| 4.30 m            | 154 kg     | 2.66 m     | 800 : 490 : 154 kg     | +7%                            |
| 4.40 m            | 146 kg     | 2.66 m     | 800 : 499 : 146 kg     | +10%                           |
| **4.50 m**        | **137 kg** | **2.28 m** | **685 : 635 : 137 kg** | **+10%**                       |
| 4.60 m and deeper | ≤ 129 kg   | none       |                        | more than 10% at any SS length |

- **EFCu** can start 0.43 m lower, at 4.50 m: 137 kg instead of 175 (−22%). Deeper, the Cu that would replace it
  sits next to the detectors at 3-7 × EFCu's activity, and the BI rises more than 10% whatever the SS.
- **SS** can then come down to 2.28 m: 685 kg instead of 619 (+66 kg). With the EFCu as built it could reach 2.69 m.
- **The SS as built** never reaches the detectors: no ROI event at any cut level in its 360 M decays.

![the designs](output/slide_designs.png)

The designs near the tube as built, coloured by their BI over the as-built's: the designs left of and below the bold
+10% line keep the BI. The star is the tube as built, the square the answer. `/Applications/LEGEND.app`'s rt shows
the same scan (`output/designs.json`) on the GDML's tube in 3D: two sliders for L1 and L2 recolour the tube by its BI
and give the masses. To build and install it:

```bash
~/venvs/v/bin/python app/make_app.py
```

## Future work

**The Hall C geometry.** Matthew Busch's string and lock update (29 Sep 2026) proposes 12-detector strings, for
running beyond 1 t, at the same central depth in the UGLAr:

- **Array**: 42 strings × 8 detectors × 3 kg = 1008 kg becomes 48 × 12 × 3 kg = 1728 kg. Each string's germanium
  grows from 1.08 to 1.64 m and its fibres from 1.36 to 1.92 m, reaching 0.28 m higher and 0.28 m deeper.
- **Tube**: 0.2 m longer, the extra height in the Cu section, and a deeper bottom head (467 → 575 mm deep): more UGLAr.
- **Lab**: the cleanroom floor, cryostat and water tank drop 0.2 m; the lock pipe grows from 7.54 to 7.81 m. It
  still waits on the crane and cleanroom details.

He makes the change conditional on the background simulations, which is what this study can answer. The array
reaches 0.28 m higher, so the EFCu likely has to start higher too; the design scan says how much. It needs
the Hall C geometry from legend-pygeom-l1000 (`units: n: 12` in its `configs/string.yaml`, 48 strings, the longer
tube and head), then the same NERSC runs and analysis on it.
