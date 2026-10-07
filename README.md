# RT backgrounds

Which mix of stainless steel (SS), OFHC copper (Cu) and electroformed copper (EFCu) the re-entrant tube
(RT) can be built from while its background stays under budget. The tube is the three materials from
the top: steel down to the seam `L1`, Cu down to `L2`, EFCu below. EFCu is the scarce one, so the
search is: **the least EFCu (`L2` as deep as possible), then the most steel (`L1` as deep as
possible)**, under a budget of `1e-5 cts/(keV·kg·yr)` at `Qbb = 2039 keV`. The answer is the
SS : Cu : EFCu mass split of that design.

**What is simulated.** remage loads the whole LEGEND-1000 geometry, so Geant4 builds every volume with
its material: the tube, the underground argon (UGLAr) inside it, the 336 HPGe detectors, their
holders, PEN plates, fibres and cryostat. One nuclide decays at a time, uniformly in the tube wall:
**Tl-208** for the ²³²Th chain and **Bi-214** for the ²³⁸U chain, the only chain members whose gammas
reach the 2 MeV region around `Qbb`. Each decay runs in full (betas, X-rays, every gamma), Geant4 carries
the gammas through the UGLAr and everything else, and the energy deposited in each HPGe detector and in
the UGLAr is recorded. Basic physics only: no optical map, no detector response.

**The analysis.** `ana/background.C` sums each detector's deposited energy per decay and applies two
cuts: **M1**, one detector fired, and an **argon veto**, rejecting a decay that leaves more than 50 keV
in the UGLAr (about 4 photoelectrons [[3]](#references)). It counts the 360 keV background window of
MAJORANA [[1]](#references) and weights every decay by the activity of the material a design puts at
its depth, so every `(L1, L2)` comes from the same run. A design passes when the simulation's 90% upper
bound on its background is under budget. A slab with no hits counts as 2.3 hits, so a design can fail
for lack of statistics, not only for too much background.

The geometry is `KSendcap_l1kGeometry.gdml` (KS) at the repository root. Depth is measured down from
the top of the tube: 6.246 m long, detectors at 4.9-5.9 m. As built: steel to 2.062 m, OFHC Cu to
4.067 m, EFCu below.

**Contents.** [Files](#files) · [Activities](#activities) · [Run](#run) · [Results](#results) ·
[Limits](#limits) · [NERSC](#nersc) · [Future work](#future-work) · [References](#references)

## Files

| source | |
| ------ | --- |
| `sim/run.mac` | the simulation: loads the GDML, records the germanium and UGLAr steps and every decay's vertex, decays one nuclide uniformly in the wall |
| `ana/background.C` | everything after it: reads the tube from the GDML and checks it, reads and checks the runs, the background as built and every design |

`sim/run.mac` takes `-s` aliases: `GDML` the geometry; `Z`, `A` the nuclide (Tl-208: 81 208, Bi-214:
83 214); `NEV` decays; `SEED` the random stream; `VOLS` a regex over the wall sections:
`reentrancetube` (EFCu: endcap, lower wall, lid), `ofhc_cu`, `ss_316l`. `ana/background.C` takes the
runs as `file=isotope` (wildcards merge files of one isotope) and, optionally, the GDML. Its printout:

| | |
| --- | --- |
| `[3]` | the tube, as read from the GDML: seams, endcap, wall thickness, and against the unmodified `l1000.gdml` |
| `[6]` | per chain: where the decays landed (all must be in the wall), hits through the cuts, and per section the volume, mass, activity, simulated decays and their density |
| `[7]` | per chain: hits by depth, and the decays each section needs to bound steel over it at 10% of the budget |
| `[8]` | the background as built, per chain and section, through the cuts; a section with no hit gets its 90% limit |
| `[9]` | for each `L2`, the deepest `L1` that passes, with its masses and SS : Cu : EFCu split |

| output (in `output/`) | made by | holds |
| --------------------- | ------- | ----- |
| `tl208.root`, `tl208_2.root`, `bi214.root`, `bi214_2.root` | `sim/run.mac` | `stp/`: 336 trees `V<string><position>` and `undergroundlar`, one row per step (`evtid`, `particle`, `edep_in_keV`, `time_in_ns`, `xloc_in_m`, `yloc_in_m`, `zloc_in_m`; the HPGe also `dist_to_surf_in_m`); `vtx`, one row per decay (`evtid`, `time_in_ns`, `xloc_in_m`, `yloc_in_m`, `zloc_in_m`, `n_part`); `detector_origins` |
| `tube.png` | `ana/background.C` | the tube, its endcap, and the endcap/barrel junction |
| `Tl208_Bi214_background.png` | `ana/background.C` | the spectrum as built, and every design's 90% bound |

git keeps only the PNGs.

## Activities

Chain activities, from Ralph's MaterialMix slides [[3]](#references):

| µBq/kg | ²³²Th chain | ²³⁸U chain | source, as quoted there |
| ------ | ----------- | ---------- | ----------------------- |
| steel | 1000 | 2500 | Bernhard |
| OFHC Cu | 1.1 | 1.3 | the MAJORANA assay paper [[4]](#references) |
| EFCu | 0.37 | 0.19 | M. Green, CD-1 |

- **Simulated nuclides:** Tl-208 decays in 35.94% of ²³²Th chain decays (the Bi-212 branch), and Bi-214
  in every ²³⁸U chain decay, assuming each chain is in equilibrium. The weights apply these factors.
- **Alternatives:** the slides also list OFHC Cu from the CD-1 at 83 µBq/kg (²³²Th) and 1.2 mBq/kg
  (²³⁸U), and ask whether those are too high. Edgar's `survival_BI.py` uses 0.19 ± 0.10 (²³⁸U) and
  < 0.077 µBq/kg (²³²Th) for EFCu.
- **Where to change them:** `cfg::chain` in `ana/background.C`. A rerun of `ana/background.C` is enough;
  the simulation does not depend on them.

## Run

From `l1000_sim/rt/`, with remage 1.1.0 (Geant4 11.3.2) and ROOT 6.40:

```bash
export G=../../KSendcap_l1kGeometry.gdml
remage -q --ignore-warnings -t 8 -w -o output/tl208.root -s GDML=$G -s Z=81 -s A=208 -s NEV=1000000 -s SEED=1 -s VOLS='reentrancetube|ofhc_cu|ss_316l' -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/tl208_2.root -s GDML=$G -s Z=81 -s A=208 -s NEV=2000000 -s SEED=3 -s VOLS='reentrancetube|ofhc_cu|ss_316l' -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/bi214.root -s GDML=$G -s Z=83 -s A=214 -s NEV=1000000 -s SEED=2 -s VOLS='reentrancetube|ofhc_cu|ss_316l' -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/bi214_2.root -s GDML=$G -s Z=83 -s A=214 -s NEV=2000000 -s SEED=4 -s VOLS='reentrancetube|ofhc_cu|ss_316l' -- sim/run.mac
root -l -b -q 'ana/background.C("output/tl208*.root=Tl208,output/bi214*.root=Bi214")'
```

| step | time (8 threads) | file |
| ---- | ---------------- | ---- |
| Tl-208, 1 M decays, seed 1 | 802 s | `tl208.root`, 32 MB |
| Tl-208, 2 M decays, seed 3 | 1526 s | `tl208_2.root`, 60 MB |
| Bi-214, 1 M decays, seed 2 | 673 s | `bi214.root`, 26 MB |
| Bi-214, 2 M decays, seed 4 | ~1190 s (1644 s on the clock, 8 min of it the laptop asleep) | `bi214_2.root`, 50 MB |
| `ana/background.C` | 38 s | `tube.png`, `Tl208_Bi214_background.png` |

- The wildcard in the last line merges every file of a chain: 3 M decays each here. More statistics:
  more runs with their own seed and file name. Never reuse a seed: two runs on one seed decay at the
  same points.
- `VOLS=ss_316l` (or `ofhc_cu`, or `reentrancetube`) puts every decay in one section. Such runs merge
  with whole-tube ones, because the weights use each section's measured decay density.

## Results

Run on 7 Oct 2026, exactly as above: 3 M decays per chain.

| | |
| --- | --- |
| background as built, after M1 + argon | **3.8 ± 1.6 × 10⁻⁸ cts/(keV·kg·yr), 0.0038 × goal**: 6 window events, all in the EFCu |
| after M1 only | 4.24 × 10⁻⁶, 0.42 × goal |
| the argon veto | keeps 1 in 44 M1 hits for Tl-208, 1 in 24 for Bi-214 |
| KS: SS : Cu : EFCu mass | 661 : 673 : 170 kg = 44 : 45 : 11 % |
| KS as a design, MC 90% bound | 24.8 × budget |
| least EFCu, then most steel | none can be shown to pass yet |
| decays the steel section needs | 4.4 × 10⁷ (Tl-208), 3.1 × 10⁸ (Bi-214) |

No decay in the steel section (0-2.06 m) reached a detector: no hit came from shallower than 2.5 m.
So every design with steel counts its steel slab at the 90% bound of a slab with no hit, mostly
Bi-214's (2.4 × 10⁻⁴, 24 × the budget), until that section has the decays above ([NERSC](#nersc)).
The background measured as built is far below the goal; the steel is unresolved, not shown to be bad.

```
[3] the tube in ../../KSendcap_l1kGeometry.gdml: z -1259.0 .. 4987.0 mm (6.246 m), 315 points, seams EFCu | 920.0 | OFHC | 2925.0 | SS
    clean surface: on-axis points 2 (expect 2), duplicates 0, spikes 0, self-intersections 0 (expect 0)
    endcap: radius 999.9998 at z -607.1, barrel 999.0000: step +0.9998 mm  <- NOT flush
    shells against the barrel: OFHC -0.0000, SS -0.0000 mm
    wall thickness [mm]:  EFCu 1.50 (z -607)  EFCu 1.50 (z 156)  EFCu 3.04 (z 910)  OFHC 6.00 (z 1922)  SS 6.00 (z 3956)
    against the unmodified ../../l1000.gdml: bottom +68.0, top +68.0, length +0.0, seams +68.0 / +68.0, barrel radius +68.0, endcap radius +69.0 mm
    tube: FAIL

window   : 1950-2350 keV minus 10 keV around 2039, 2103.5, 2118.5, 2204.1 keV (MAJORANA's BEW): 360 keV
cuts     : M1 (exactly one detector above 5 keV), argon veto (UGLAr deposit <= 50 keV). no detector model: deposited energy
activity : chain [uBq/kg] 232Th / 238U: steel 1000 / 2500, Cu 1.1 / 1.3, EFCu 0.37 / 0.19; Tl208 is 35.94% of 232Th, Bi214 100% of 238U

=== Tl208   (output/tl208.root+output/tl208_2.root, 3000000 decays)
[6] decays: EFCu 18.0%, OFHC 40.1%, SS 41.9%; in the argon 0, outside the wall 0   (merged from 2 files)
    hits: 30289, M1 25213, M1 + argon 574; in the window: 968, 785, 2
    section mat      V [m^3]   M [kg]     A [Bq]   MC dec.  MC per m^3  decays/yr per
    mother  EFCu     0.02535    226.4  3.011e-05    539165    21266202        0.00176
    OFHC    Cu       0.07496    671.7  2.655e-04   1203540    16055198        0.00696
    SS      steel    0.07742    611.7  2.198e-01   1257295    16239009           5.52
    sampling density mother / shells = 1.317   (1.000 would be uniform; the weights use the measured density)
[7] depth [m]          decays     hits    hits/decay   window
     0.00..0.62         513177        0      0.00e+00        0
     0.62..1.25         380435        0      0.00e+00        0
     1.25..1.87         380637        0      0.00e+00        0
     1.87..2.50         377795        0      0.00e+00        0
     2.50..3.12         377352        5      1.33e-05        0
     3.12..3.75         375453       26      6.92e-05        2
     3.75..4.37         256402      442      1.72e-03       12
     4.37..5.00         124889     6138      4.91e-02      199
     5.00..5.62         125378    14719      1.17e-01      452
     5.62..6.25          88482     8959      1.01e-01      303   <- nearest the detectors
    a slab with no window hit is bounded at 2.30 x one decay's weight. to bound steel there at 10% of the budget:
      steel over SS        1257295 decays now: bound 3.5e-05, needs 4.4e+07 decays in it
      steel over OFHC      1203540 decays now: bound 3.6e-05, needs 4.3e+07 decays in it
      steel over mother     539165 decays now: bound 2.7e-05, needs 1.5e+07 decays in it

=== Bi214   (output/bi214.root+output/bi214_2.root, 3000000 decays)
[6] decays: EFCu 17.9%, OFHC 40.1%, SS 41.9%; in the argon 0, outside the wall 0   (merged from 2 files)
    hits: 15587, M1 13810, M1 + argon 582; in the window: 50, 48, 4
    section mat      V [m^3]   M [kg]     A [Bq]   MC dec.  MC per m^3  decays/yr per
    mother  EFCu     0.02535    226.4  4.302e-05    538432    21237291        0.00252
    OFHC    Cu       0.07496    671.7  8.732e-04   1203167    16050222         0.0229
    SS      steel    0.07742    611.7  1.529e+00   1258401    16253294           38.3
    sampling density mother / shells = 1.315   (1.000 would be uniform; the weights use the measured density)
[7] depth [m]          decays     hits    hits/decay   window
     0.00..0.62         512700        0      0.00e+00        0
     0.62..1.25         381669        0      0.00e+00        0
     1.25..1.87         380617        0      0.00e+00        0
     1.87..2.50         376712        0      0.00e+00        0
     2.50..3.12         377397        0      0.00e+00        0
     3.12..3.75         376533        9      2.39e-05        0
     3.75..4.37         256453      136      5.30e-04        0
     4.37..5.00         124415     3059      2.46e-02       12
     5.00..5.62         124831     7552      6.05e-02       25
     5.62..6.25          88673     4831      5.45e-02       13   <- nearest the detectors
    a slab with no window hit is bounded at 2.30 x one decay's weight. to bound steel there at 10% of the budget:
      steel over SS        1258401 decays now: bound 2.4e-04, needs 3.1e+08 decays in it
      steel over OFHC      1203167 decays now: bound 2.5e-04, needs 3.0e+08 decays in it
      steel over mother     538432 decays now: bound 1.9e-04, needs 1.0e+08 decays in it

[8] background index as built [cts/(keV kg yr)], +- MC statistics; a section with no window hit gets its 90% limit
    chain  sect       win hits   no cut                 M1                     M1 + argon            
    Tl208  steel         0/0/0   < 3.45e-05 (90%)       < 3.45e-05 (90%)       < 3.45e-05 (90%)      
           Cu            5/4/0   9.67e-08 +- 4.3e-08    7.74e-08 +- 3.9e-08    < 4.44e-08 (90%)      
           EFCu      963/781/2   4.71e-06 +- 1.5e-07    3.82e-06 +- 1.4e-07    9.79e-09 +- 6.9e-09   
           tube                  4.81e-06 +- 1.6e-07    3.90e-06 +- 1.4e-07    9.79e-09 +- 6.9e-09   
    Bi214  steel         0/0/0   < 2.40e-04 (90%)       < 2.40e-04 (90%)       < 2.40e-04 (90%)      
           Cu            0/0/0   < 1.46e-07 (90%)       < 1.46e-07 (90%)       < 1.46e-07 (90%)      
           EFCu        50/48/4   3.50e-07 +- 5.0e-08    3.36e-07 +- 4.9e-08    2.80e-08 +- 1.4e-08   
           tube                  3.50e-07 +- 5.0e-08    3.36e-07 +- 4.9e-08    2.80e-08 +- 1.4e-08   
    ALL CHAINS   no cut: 5.16e-06 +- 1.7e-07 (0.52 x goal)   M1: 4.24e-06 +- 1.5e-07 (0.42 x goal)   M1 + argon: 3.78e-08 +- 1.6e-08 (0.0038 x goal)
    win hits: no cut / M1 / M1 + argon. sections with no window hit are left out of the totals

[9] designs: steel to L1, Cu to L2, EFCu below. MC alone, after M1 + argon; passes if its 90% bound <= the budget 1e-05
    L2 [m]  EFCu [kg]  steel to steel [kg]   Cu [kg]   SS:Cu:EFCu mass          BI   90% bound  x budget
    3.00          526         -   no steel can be shown to pass yet
    3.25          442         -   no steel can be shown to pass yet
    3.50          358         -   no steel can be shown to pass yet
    3.75          275         -   no steel can be shown to pass yet
    4.00          191         -   no steel can be shown to pass yet
    4.07          170    2.06 m        661       673    44 : 45 : 11 %    3.78e-08    2.48e-04     24.83   <- KS as built
    4.25          152         -   no steel can be shown to pass yet
    4.50          131         -   no steel can be shown to pass yet
    4.75          110         -   no steel can be shown to pass yet
    5.00           89         -   no steel can be shown to pass yet
    5.25           68         -   no steel can be shown to pass yet
    5.50           47         -   no steel can be shown to pass yet
    5.75           26         -   no steel can be shown to pass yet
    6.00            8         -   no steel can be shown to pass yet
    no design with steel passes on the MC alone: [7] says how many decays its slabs need
    masses in kg; steel includes the lid, which sits in the top slab

wrote output/Tl208_Bi214_background.png and output/tube.png
```

![background](output/Tl208_Bi214_background.png)

Left: the spectrum as built, through the cuts; the argon veto removes almost everything, and the
2614.5 keV line of Tl-208 is what survives longest. Right: every design's 90% bound in units of the
budget; the star is KS. The flat ~25 × wherever there is steel is the bound of the empty steel slab,
not a measured background.

![the tube](output/tube.png)

`[3]` says `FAIL` for one reason: the KS endcap is 1.0 mm wider than the barrel where they meet, at
z = −607 mm (right panel). The walls: EFCu 1.5 mm, thickening to 3 mm near the OFHC seam; OFHC and
steel 6 mm. Against the unmodified geometry (dashed), everything moved up 68 mm and out 68 mm.

`[6]` checks the decays: all are in the wall. remage 1.1 fills the EFCu (the mother volume) about 1.3 ×
as densely as the two shells inside it; the weights use each section's measured density, so this
costs nothing but statistics.

## Limits

- No detector model: the energy is what Geant4 deposited.
- The argon veto counts deposited energy, as if every deposit above 50 keV were seen: real light
  collection (the optical map) would veto less.
- Reweighting changes activities, not geometry: it holds while the tube and the detectors stay put.
- The lid at the top takes the top slab's material in every design, KS included.
- Each chain is assumed in equilibrium down to Tl-208 and Bi-214.
- A design with steel can only pass once its steel slab has enough decays.

## NERSC

The decays `[7]` asks for, as job arrays of 10⁷: stage 1 fills the steel section, stage 2 the OFHC
section, where steel would reach below the seam. The task counts come from the run above; resize them
to your `[7]`. Once, from here, with your login and path:

```bash
ssh <user>@perlmutter.nersc.gov "mkdir -p <path>/LEGEND1000-Simulation/l1000_sim/rt/output <path>/LEGEND1000-Simulation/l1000_sim/rt/sim"
rsync -av ../../KSendcap_l1kGeometry.gdml <user>@perlmutter.nersc.gov:<path>/LEGEND1000-Simulation/
rsync -av sim/ <user>@perlmutter.nersc.gov:<path>/LEGEND1000-Simulation/l1000_sim/rt/sim/
```

On Perlmutter, from `l1000_sim/rt/` (`m2676` is assumed to be the LEGEND allocation):

```bash
shifterimg pull docker:legendexp/remage:v1.1.0
# stage 1, the steel section
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --output=output/%x_%a.log --image=docker:legendexp/remage:v1.1.0 --array=3001-3005 --job-name=rt_tl208_ss --wrap "shifter remage -q --ignore-warnings -t 32 -w -s GDML=../../KSendcap_l1kGeometry.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/tl208_\$SLURM_ARRAY_TASK_ID.root -s Z=81 -s A=208 -s VOLS=ss_316l -- sim/run.mac"
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --output=output/%x_%a.log --image=docker:legendexp/remage:v1.1.0 --array=4001-4031 --job-name=rt_bi214_ss --wrap "shifter remage -q --ignore-warnings -t 32 -w -s GDML=../../KSendcap_l1kGeometry.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/bi214_\$SLURM_ARRAY_TASK_ID.root -s Z=83 -s A=214 -s VOLS=ss_316l -- sim/run.mac"
# stage 2, the OFHC section
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --output=output/%x_%a.log --image=docker:legendexp/remage:v1.1.0 --array=5001-5005 --job-name=rt_tl208_cu --wrap "shifter remage -q --ignore-warnings -t 32 -w -s GDML=../../KSendcap_l1kGeometry.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/tl208_\$SLURM_ARRAY_TASK_ID.root -s Z=81 -s A=208 -s VOLS=ofhc_cu -- sim/run.mac"
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --output=output/%x_%a.log --image=docker:legendexp/remage:v1.1.0 --array=6001-6030 --job-name=rt_bi214_cu --wrap "shifter remage -q --ignore-warnings -t 32 -w -s GDML=../../KSendcap_l1kGeometry.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/bi214_\$SLURM_ARRAY_TASK_ID.root -s Z=83 -s A=214 -s VOLS=ofhc_cu -- sim/run.mac"
squeue --me
```

Each task takes about 40 minutes and writes ~300 MB. After each stage, here:

```bash
rsync -av --include='*_[0-9]*.root' --exclude='*' <user>@perlmutter.nersc.gov:<path>/LEGEND1000-Simulation/l1000_sim/rt/output/ output/
root -l -b -q 'ana/background.C("output/tl208*.root=Tl208,output/bi214*.root=Bi214")'
```

## Future work

| ‹placeholder› | adds | needs |
| ------------- | ---- | ----- |
| ‹optical map› | a realistic argon veto: the light the SiPMs see, not the energy deposited | the LAr optical map |
| ‹detector response› | dead layer and resolution | each detector's FCCD and resolution |
| ‹PSD› | pulse-shape discrimination (A/E) | drift-time maps for each crystal |
| ‹survival check› | cut survival against the design report and Edgar's chain | the cuts above |
| ‹reach extrapolation› | an estimate where the MC has no hits, without more decays | a validated model of how the hit probability falls with height |
| ‹per-detector rates› | the tube's fingerprint in data: each detector's 2615 and 1764 keV rates | nothing new |
| ‹purity requirement› | the activity at which a section alone gives the goal | nothing new |
| ‹Hall C geometry› | 12-detector strings and a longer tube | a new GDML |
| ‹margin› | a design that survives activities being off (MAJORANA's came out ~5 x low [[1]](#references)) | a budget below the whole goal |

## References

1. C.R. Haufe et al. (MAJORANA), "Modeling Backgrounds for the MAJORANA DEMONSTRATOR", arXiv:2209.10592 (2023).
2. I.J. Arnquist et al. (MAJORANA), "Final result of the MAJORANA DEMONSTRATOR's search for neutrinoless
   double-β decay in ⁷⁶Ge", arXiv:2207.07638 (2022): the window's three excluded lines.
3. R. Massarczyk, "Re-entrant tube material combination vs ROI", 23 Jun 2026: `l1000_sim/Ralphs_sim/2026-06-23-MaterialMix.pdf`
   (activities on pp. 5 and 15, the argon threshold on p. 4).
4. MAJORANA Collaboration, "Assay-based background projection for the MAJORANA DEMONSTRATOR using Monte
   Carlo uncertainty propagation", Phys. Rev. C 110, 055804 (2024).
