# RT backgrounds

The cheapest re-entrant tube (RT) whose background stays under budget. The tube is three materials
from the top: steel down to the seam `L1`, OFHC Cu down to `L2`, electroformed copper (EFCu) below.
EFCu is the cleanest and the most expensive, so the search is: **the least EFCu (`L2` as deep as
possible), then the most steel (`L1` as deep as possible)**, under a budget of `1e-5 cts/(keV·kg·yr)`
at `Qbb = 2039 keV`.

**Method.** remage decays Tl-208 (²³²Th chain) or Bi-214 (²³⁸U chain) uniformly through the tube wall
and records every step in the germanium. `ana/background.C` sums each detector's deposited energy,
keeps the decays where one detector fired (M1), and counts the 360 keV background window of MAJORANA
[[1]](#references). Each decay is weighted by the activity of the material a design puts at its
depth, so every `(L1, L2)` comes from the same run. A design passes when the MC's 90% upper bound is
under budget.

The geometry is `KSendcap_l1kGeometry.gdml` (KS) at the repository root. Depth is measured down from
the top of the tube: 6.246 m long, detectors at 4.9-5.9 m. As built: steel to 2.062 m, OFHC Cu to
4.067 m, EFCu below.

**Contents.** [Files](#files) · [Run](#run) · [Results](#results) · [Limits](#limits) ·
[NERSC](#nersc) · [Future work](#future-work) · [References](#references)

## Files

| source | |
| ------ | --- |
| `sim/source.mac` | one nuclide (`Z`, `A`) at rest, uniformly in the wall sections `VOLS`, `NEV` decays, seed `SEED` |
| `sim/run.mac` | loads the GDML, records the germanium steps and every decay's vertex, runs `source.mac` |
| `sim/rt.h` | reads the tube's shape from the GDML: its sections, volumes and the material at any point |
| `sim/tube.C` | is the tube one clean surface, and how KS differs from the unmodified `l1000.gdml` |
| `ana/vertices.C` | did the decays fill the wall, and only the wall |
| `ana/background.C` | the background as built, and every steel / Cu / EFCu design |

| output (in `output/`) | made by | holds |
| --------------------- | ------- | ----- |
| `tl208.root`, `tl208_2.root`, `bi214.root`, `bi214_2.root` | remage, `sim/run.mac` | `stp/`: 336 trees `V<string><position>`, one row per germanium step (`evtid`, `particle`, `edep_in_keV`, `time_in_ns`, `xloc_in_m`, `yloc_in_m`, `zloc_in_m`, `dist_to_surf_in_m`); `vtx`, one row per decay (`evtid`, `time_in_ns`, `xloc_in_m`, `yloc_in_m`, `zloc_in_m`, `n_part`); `detector_origins` |
| `tube.png` | `sim/tube.C` | the tube, its endcap, and the endcap/barrel junction |
| `Tl208_Bi214_background.png` | `ana/background.C` | the spectrum as built, and every design's 90% bound |

git keeps only the PNGs.

## Run

From `l1000_sim/rt/`, with remage v0.26.0 and ROOT 6.40:

```bash
export G=../../KSendcap_l1kGeometry.gdml
root -l -b -q sim/tube.C
remage -q --ignore-warnings -t 8 -w -o output/tl208.root -s GDML=$G -s Z=81 -s A=208 -s NEV=1000000 -s SEED=1 -s VOLS='reentrancetube|ofhc_cu|ss_316l' -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/tl208_2.root -s GDML=$G -s Z=81 -s A=208 -s NEV=2000000 -s SEED=3 -s VOLS='reentrancetube|ofhc_cu|ss_316l' -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/bi214.root -s GDML=$G -s Z=83 -s A=214 -s NEV=1000000 -s SEED=2 -s VOLS='reentrancetube|ofhc_cu|ss_316l' -- sim/run.mac
remage -q --ignore-warnings -t 8 -w -o output/bi214_2.root -s GDML=$G -s Z=83 -s A=214 -s NEV=2000000 -s SEED=4 -s VOLS='reentrancetube|ofhc_cu|ss_316l' -- sim/run.mac
root -l -b -q 'ana/vertices.C("output/tl208.root")'
root -l -b -q 'ana/background.C("output/tl208*.root=Tl208,output/bi214*.root=Bi214")'
```

| step | time (8 threads) | file |
| ---- | ---------------- | ---- |
| `sim/tube.C` | < 1 s | `tube.png` |
| Tl-208, 1 M decays, seed 1 | 633 s | `tl208.root`, 25 MB |
| Tl-208, 2 M decays, seed 3 | 1222 s | `tl208_2.root`, 48 MB |
| Bi-214, 1 M decays, seed 2 | 570 s | `bi214.root`, 23 MB |
| Bi-214, 2 M decays, seed 4 | 1098 s | `bi214_2.root`, 45 MB |
| `ana/vertices.C` | 6 s | |
| `ana/background.C` | 26 s | `Tl208_Bi214_background.png` |

- `sim/tube.C` exits with code 1: it reports `FAIL`, as expected (see below).
- The wildcard in the last line merges every file of a chain: 3 M decays each here. More statistics:
  more runs with their own seed and file name. Never reuse a seed: two runs on one seed decay at the
  same points.
- `VOLS=ss_316l` (or `ofhc_cu`, or `reentrancetube`, the EFCu) puts every decay in one section. Such
  runs merge with whole-tube ones, because the weights use each section's measured decay density.

## Results

Run on 6 Oct 2026, exactly as above: 3 M decays per chain.

| | |
| --- | --- |
| background as built, after M1 | **< 2.59 × 10⁻⁶ cts/(keV·kg·yr), < 0.26 × goal** (< 3.02 × 10⁻⁶ without M1) |
| where it comes from | Tl-208 in the EFCu: 501 of the 535 window hits after M1 |
| KS as a design, MC 90% bound | 28.8 × budget |
| least EFCu, then most steel | none can be shown to pass yet |
| decays the steel section needs | 1.2 × 10⁸ (Tl-208), 3.1 × 10⁸ (Bi-214) |

The `<` is there because EFCu's Tl-208 activity is an upper limit. No decay in the steel section
(0-2.06 m) reached a detector: no hit came from shallower than 2.5 m. So every design with steel
counts its steel slab at the 90% bound of a slab with no hit, about 30 × the budget, until that
section has the decays above ([NERSC](#nersc)).

### The tube — `sim/tube.C`

```
file : ../../KSendcap_l1kGeometry.gdml
RT   : z -1259.0 .. 4987.0 mm, 315 points, seams EFCu | 920.0 | OFHC | 2925.0 | SS

[1] solids:
      reentrancetube         n=315  r    0.0000 ..  999.9998   z  -1259.00 ..   4987.00
      undergroundlar         n=315  r    0.0000 ..  998.4983   z  -1258.99 ..   4987.00
      ofhc_cu_outer_bound    n=44   r    0.0000 ..  999.0000   z    920.00 ..   2925.00
      ofhc_cu_inner_bound    n=44   r    0.0000 ..  997.5000   z    920.00 ..   2925.00
      ss_316l_outer_bound    n=45   r    0.0000 ..  999.0000   z   2925.00 ..   4987.00
      ss_316l_inner_bound    n=45   r    0.0000 ..  993.0000   z   2925.00 ..   4987.00

[2] clean surface: on-axis points 2 (expect 2), duplicates 0, spikes 0, self-intersections 0 (expect 0)

[3] endcap flush: barrel r 999.0000 (z 1922.5), endcap r 999.9998 (z -607.1), step +0.9998 mm (tolerance 0.050)  <- NOT flush

[4] shells:
      ofhc_cu_outer_bound    rmax  999.0000 (barrel -0.0000)
      ss_316l_outer_bound    rmax  999.0000 (barrel -0.0000)

[5] wall thickness:
          z [mm]   sect      r_out       r_in   t [mm]
          -607.1   EFCu   999.9998   998.4983   1.5015
           156.4   EFCu   999.0000   997.5000   1.5000
           910.0   EFCu   999.0000   995.9580   3.0420
          1922.5   OFHC   999.0000   993.0000   6.0000
          2915.0   OFHC   999.0000   993.0000   6.0000
          3956.0     SS   999.0000   993.0000   6.0000
          4977.0     SS   999.0000   993.0000   6.0000

[6] KS vs mint ../../l1000.gdml:
                              mint          KS      delta
      RT z bottom       -1327.0000  -1259.0000   +68.0000
      RT z top           4919.0000   4987.0000   +68.0000
      RT length          6246.0000   6246.0000    +0.0000
      endcap rmax         931.0000    999.9998   +68.9998
      barrel r            931.0000    999.0000   +68.0000
      EFCu|OFHC seam      852.0000    920.0000   +68.0000
      OFHC|SS seam       2857.0000   2925.0000   +68.0000

RESULT [tube]: FAIL
```

![the tube](output/tube.png)

`FAIL` for one reason: `[3]`, the KS endcap is 1.0 mm wider than the barrel where it meets it, at
z = −607 mm (right panel). The walls: EFCu 1.5 mm, thickening to 3 mm near the OFHC seam; OFHC and steel
6 mm. Against the unmodified geometry (dashed), everything moved up 68 mm and out 68 mm.

### The decays — `ana/vertices.C`

```
file     : output/tl208.root
geometry : ../../KSendcap_l1kGeometry.gdml
vertices : 1000000    z -1255.2 .. 4987.0 mm    r 3.5 .. 1000.0 mm
RT extent: -1259.0 .. 4987.0 mm    seams 920.0 / 2925.0

[2] material (exact point-in-solid, not a z-cut):
      EFCu   117754 ( 11.8%)
      OFHC   433549 ( 43.4%)
      SS     448697 ( 44.9%)

[3] coverage along z:
      -1259..    -869 mm | ###                                      9314
       -869..    -478 mm | #####                                    16933
       -478..     -88 mm | #####                                    17094
        -88..     302 mm | ######                                   17144
        302..     693 mm | #####                                    16893
        693..    1083 mm | ###############                          45290
       1083..    1474 mm | #############################            84589
       1474..    1864 mm | #############################            84621
       1864..    2254 mm | #############################            85091
       2254..    2645 mm | #############################            84755
       2645..    3035 mm | #############################            85039
       3035..    3426 mm | #############################            84953
       3426..    3816 mm | #############################            84631
       3816..    4206 mm | #############################            84554
       4206..    4597 mm | #############################            85078
       4597..    4987 mm | ######################################## 114021

RESULT [vertices]: PASS
```

Every vertex is in the wall. A fill uniform by volume would put 14.3% in the EFCu; remage fills the
EFCu mother about 20% more sparsely than the shells inside it, which `ana/background.C` corrects with
each section's measured decay density (`[5]`). The fuller last slab is the lid; the thinner slabs
below 693 mm are the 1.5 mm EFCu wall.

### The background — `ana/background.C`

```
geometry : ../../KSendcap_l1kGeometry.gdml
RT wall  : 0.17774 m^3 over 6.246 m
window   : 1950-2350 keV minus 10 keV around 2039, 2103.5, 2118.5, 2204.1 keV (MAJORANA's BEW): 360 keV
cut      : M1 (exactly one detector above 5 keV). no detector model: deposited energy

=== Tl208   (output/tl208.root+output/tl208_2.root, 3000000 decays)
[4] hits: 19836, after M1 16717; in the window: 588, after M1 503   (merged from 2 files)
[5] normalisation:  volume  mat      V [m^3]   M [kg]     A [Bq]   MC dec.  MC per m^3  decays/yr per
                    mother  EFCu     0.02535    226.4  1.743e-05    353864    13957403        0.00155
                    OFHC    Cu       0.07496    671.7  6.717e-04   1300295    17345908         0.0163
                    SS      steel    0.07742    611.7  6.117e-01   1345841    17382654           14.3
    sampling density mother / shells = 0.804   (1.000 would be uniform; the weights use the measured density)
[6] depth [m]          decays     hits    hits/decay   window
     0.00..0.62         494859        0      0.00e+00        0
     0.62..1.25         406996        0      0.00e+00        0
     1.25..1.87         407429        0      0.00e+00        0
     1.87..2.50         407429        0      0.00e+00        0
     2.50..3.12         407569        2      4.91e-06        0
     3.12..3.75         405600       34      8.38e-05        0
     3.75..4.37         247869      334      1.35e-03       11
     4.37..5.00          81988     4053      4.94e-02      120
     5.00..5.62          82171     9532      1.16e-01      277
     5.62..6.25          58090     5881      1.01e-01      180   <- nearest the detectors
    a slab with no window hit is bounded at 2.30 x one decay's weight. to bound steel there at 10% of the budget:
      steel over SS        1345841 decays now: bound 9.2e-05, needs 1.2e+08 decays in it
      steel over OFHC      1300295 decays now: bound 9.2e-05, needs 1.2e+08 decays in it
      steel over mother     353864 decays now: bound 1.1e-04, needs 4.0e+07 decays in it

=== Bi214   (output/bi214.root+output/bi214_2.root, 3000000 decays)
[4] hits: 10290, after M1 9023; in the window: 34, after M1 32   (merged from 2 files)
[5] normalisation:  volume  mat      V [m^3]   M [kg]     A [Bq]   MC dec.  MC per m^3  decays/yr per
                    mother  EFCu     0.02535    226.4  4.302e-05    354125    13967698        0.00383
                    OFHC    Cu       0.07496    671.7  6.717e-04   1299744    17338558         0.0163
                    SS      steel    0.07742    611.7  1.529e+00   1346131    17386400           35.8
    sampling density mother / shells = 0.804   (1.000 would be uniform; the weights use the measured density)
[6] depth [m]          decays     hits    hits/decay   window
     0.00..0.62         493509        0      0.00e+00        0
     0.62..1.25         408551        0      0.00e+00        0
     1.25..1.87         407184        0      0.00e+00        0
     1.87..2.50         406029        0      0.00e+00        0
     2.50..3.12         407424        0      0.00e+00        0
     3.12..3.75         406678        5      1.23e-05        0
     3.75..4.37         248094       91      3.67e-04        1
     4.37..5.00          81900     2037      2.49e-02        6
     5.00..5.62          82083     4964      6.05e-02       16
     5.62..6.25          58548     3193      5.45e-02       11   <- nearest the detectors
    a slab with no window hit is bounded at 2.30 x one decay's weight. to bound steel there at 10% of the budget:
      steel over SS        1346131 decays now: bound 2.3e-04, needs 3.1e+08 decays in it
      steel over OFHC      1299744 decays now: bound 2.3e-04, needs 3.0e+08 decays in it
      steel over mother     354125 decays now: bound 2.9e-04, needs 1.0e+08 decays in it

[7] background index as built [cts/(keV kg yr)]: BI +- MC statistics; a section with no window hit gets its 90% limit
    chain  sect    win hits   no cut                     M1                        
    Tl208  steel     0/0      < 9.30e-05 (90%)           < 9.30e-05 (90%)          
           Cu        2/3      1.36e-07 +- 7.8e-08        9.06e-08 +- 6.4e-08       
           EFCu    501/585    <2.53e-06 +- 1.0e-07       <2.16e-06 +- 9.7e-08        activity is an upper limit
           tube               <2.66e-06 +- 1.3e-07       <2.25e-06 +- 1.2e-07      
    Bi214  steel     0/0      < 2.32e-04 (90%)           < 2.32e-04 (90%)          
           Cu        0/0      < 1.04e-07 (90%)           < 1.04e-07 (90%)          
           EFCu     32/34     3.62e-07 +- 6.2e-08        3.41e-07 +- 6.0e-08       
           tube               3.62e-07 +- 6.2e-08        3.41e-07 +- 6.0e-08       
    ALL CHAINS          <3.02e-06 +- 1.4e-07 (<0.30 x goal)   M1: <2.59e-06 +- 1.3e-07 (<0.26 x goal)
    sections with no window hit are left out of the totals; win hits reads M1/no cut

[8] designs: steel to L1, Cu to L2, EFCu below. MC alone, after M1; passes if its 90% bound <= the budget 1e-05
    L2 [m]  EFCu [kg]  steel to steel [kg]   Cu [kg]          BI   90% bound  x budget
    3.00          526         -   no steel can be shown to pass yet
    3.25          442         -   no steel can be shown to pass yet
    3.50          358         -   no steel can be shown to pass yet
    3.75          275         -   no steel can be shown to pass yet
    4.00          191         -   no steel can be shown to pass yet
    4.07          170    2.06 m        661       673    2.59e-06    2.88e-04     28.80   <- KS as built
    4.25          152         -   no steel can be shown to pass yet
    4.50          131         -   no steel can be shown to pass yet
    4.75          110         -   no steel can be shown to pass yet
    5.00           89         -   no steel can be shown to pass yet
    5.25           68         -   no steel can be shown to pass yet
    5.50           47         -   no steel can be shown to pass yet
    5.75           26         -   no steel can be shown to pass yet
    6.00            8         -   no steel can be shown to pass yet
    no design with steel passes on the MC alone: [6] says how many decays its slabs need
    masses in kg; steel includes the lid, which sits in the top slab

wrote output/Tl208_Bi214_background.png
```

| | |
| --- | --- |
| `[4]` | hits and window hits, no cut and M1 |
| `[5]` | per section: volume, mass, activity, simulated decays, their density, real decays per year per simulated one |
| `[6]` | hits by depth, and the decays each section needs to bound steel over it at 10% of the budget |
| `[7]` | the background as built, per chain and section; a section with no window hit gets its 90% limit |
| `[8]` | for each `L2`, the deepest `L1` whose 90% bound passes, with its masses |

Activities in µBq/kg (Bi-214 / Tl-208): steel 2500 / 1000 and Cu 1 / 1 (Ralph's materialMix); EFCu
0.19 / < 0.077 (Edgar's radioassay).

![background](output/Tl208_Bi214_background.png)

Left: the spectrum as built, without and with M1; the 2614.5 keV line of Tl-208 stands out. Right:
every design's 90% bound in units of the budget; the star is KS. The flat ~30 × wherever there is
steel is the bound of the empty steel slab, not a measured background.

## Limits

- No detector model: the energy is what Geant4 deposited. The only cut is M1, so every number is an
  upper estimate of the background after all cuts.
- Reweighting changes activities, not geometry: it holds while the tube and the detectors stay put.
- The lid at the top takes the top slab's material in every design, KS included.
- A design with steel can only pass once its steel slab has enough decays.

## NERSC

The decays `[6]` asks for, as job arrays of 10⁷: stage 1 fills the steel section, stage 2 the OFHC
section, where steel would reach below the seam. The task counts come from a 3 M-decay run; resize
them to your `[6]`. Once, from here, with your login and path:

```bash
ssh <user>@perlmutter.nersc.gov "mkdir -p <path>/LEGEND1000-Simulation/l1000_sim/rt/output <path>/LEGEND1000-Simulation/l1000_sim/rt/sim"
rsync -av ../../KSendcap_l1kGeometry.gdml <user>@perlmutter.nersc.gov:<path>/LEGEND1000-Simulation/
rsync -av sim/ <user>@perlmutter.nersc.gov:<path>/LEGEND1000-Simulation/l1000_sim/rt/sim/
```

On Perlmutter, from `l1000_sim/rt/` (`m2676` is assumed to be the LEGEND allocation):

```bash
shifterimg pull docker:legendexp/remage:v0.26.0
# stage 1, the steel section
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --output=output/%x_%a.log --image=docker:legendexp/remage:v0.26.0 --array=3001-3012 --job-name=rt_tl208_ss --wrap "shifter remage -q --ignore-warnings -t 32 -w -s GDML=../../KSendcap_l1kGeometry.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/tl208_\$SLURM_ARRAY_TASK_ID.root -s Z=81 -s A=208 -s VOLS=ss_316l -- sim/run.mac"
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --output=output/%x_%a.log --image=docker:legendexp/remage:v0.26.0 --array=4001-4031 --job-name=rt_bi214_ss --wrap "shifter remage -q --ignore-warnings -t 32 -w -s GDML=../../KSendcap_l1kGeometry.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/bi214_\$SLURM_ARRAY_TASK_ID.root -s Z=83 -s A=214 -s VOLS=ss_316l -- sim/run.mac"
# stage 2, the OFHC section
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --output=output/%x_%a.log --image=docker:legendexp/remage:v0.26.0 --array=5001-5012 --job-name=rt_tl208_cu --wrap "shifter remage -q --ignore-warnings -t 32 -w -s GDML=../../KSendcap_l1kGeometry.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/tl208_\$SLURM_ARRAY_TASK_ID.root -s Z=81 -s A=208 -s VOLS=ofhc_cu -- sim/run.mac"
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --output=output/%x_%a.log --image=docker:legendexp/remage:v0.26.0 --array=6001-6030 --job-name=rt_bi214_cu --wrap "shifter remage -q --ignore-warnings -t 32 -w -s GDML=../../KSendcap_l1kGeometry.gdml -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -o output/bi214_\$SLURM_ARRAY_TASK_ID.root -s Z=83 -s A=214 -s VOLS=ofhc_cu -- sim/run.mac"
squeue --me
```

Each task takes about an hour and writes ~200 MB. After each stage, here:

```bash
rsync -av --include='*_[0-9]*.root' --exclude='*' <user>@perlmutter.nersc.gov:<path>/LEGEND1000-Simulation/l1000_sim/rt/output/ output/
root -l -b -q 'ana/background.C("output/tl208*.root=Tl208,output/bi214*.root=Bi214")'
```

## Future work

| ‹placeholder› | adds | needs |
| ------------- | ---- | ----- |
| ‹optical map› | the argon veto | the LAr optical map, and argon steps recorded in `sim/run.mac` |
| ‹detector response› | dead layer and resolution | each detector's FCCD and resolution |
| ‹PSD› | pulse-shape discrimination (A/E) | drift-time maps for each crystal |
| ‹survival check› | cut survival against the design report and Edgar's chain | the two cuts above |
| ‹reach extrapolation› | an estimate where the MC has no hits, without more decays | a validated model of how the hit probability falls with height |
| ‹per-detector rates› | the tube's fingerprint in data: each detector's 2615 and 1764 keV rates | nothing new |
| ‹purity requirement› | the activity at which a section alone gives the goal | nothing new |
| ‹Hall C geometry› | 12-detector strings and a longer tube | a new GDML |
| ‹margin› | a design that survives activities being off (MAJORANA's came out ~5 x low [[1]](#references)) | a budget below the whole goal |

## References

1. C.R. Haufe et al. (MAJORANA), "Modeling Backgrounds for the MAJORANA DEMONSTRATOR", arXiv:2209.10592 (2023).
2. I.J. Arnquist et al. (MAJORANA), "Final result of the MAJORANA DEMONSTRATOR's search for neutrinoless
   double-β decay in ⁷⁶Ge", arXiv:2207.07638 (2022): the window's three excluded lines.
