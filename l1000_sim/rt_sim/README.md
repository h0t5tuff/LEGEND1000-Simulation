# RT backgrounds

**goal.** How much background does the re-entrant tube (RT) put into the Ge detectors,
and does the choice and placement of its materials matter?
**Target.** `1e-5 cts/(keV·kg·yr)` in the ROI, `2039 ± 55 keV`.

**sim.** Simulate each decay chain once, uniformly through the whole tube wall, then give
every decay a weight for the material a design puts at its depth. Any design — which
material, where the seams sit — becomes arithmetic on one run, not a new simulation.

One run per chain: **Tl208** (Th232) and **Bi214** (U238). The background is their sum.

## Pipeline

`geom/` the geometry, `sim/` what remage runs, `ana/` what reads its output.
`geom/rt.h` reads the tube from the GDML. no macro hard-codes a number

```
⓪ geometry   geom/geo_check.mac     (clean GDML)
             geom/geo_rt.C          (clean RT)                        -> geo_rt.png
① remage     sim/sim_run.mac + sim/src_events.mac  (decays)           -> <run>.root
             ana/src_vertices.C (decays fill) ana/sim_hits.C (Ge hit) -> <run>_hits.png
② response   ana/sim_psd.py         (dead layer + A/E)                -> <run>_psd.csv
③ cuts       ana/sim_background.C   (cuts, BI)                        -> <Isos>_background.png
                                                                         <Isos>_survival.png
```

Everything is run from `l1000_sim/rt_sim/` and writes to `output/`:

```bash
export G=../../KSendcap_l1kGeometry.gdml
```

## Run

1 M decays per chain

```bash
remage -q --ignore-warnings -s GDML=$G -s SKIP=true -s NPOINTS=1000 -- geom/geo_check.mac
```

```bash
root -l -b -q geom/geo_rt.C
```

```bash
remage -q --ignore-warnings -t 8 -w -o output/tl208.root -s GDML=$G -s Z=81 -s A=208 -s NEV=1000000 -s SEED=1 -- sim/sim_run.mac
```

```bash
remage -q --ignore-warnings -t 8 -w -o output/bi214.root -s GDML=$G -s Z=83 -s A=214 -s NEV=1000000 -s SEED=2 -- sim/sim_run.mac
```

```bash
root -l -b -q 'ana/src_vertices.C("output/tl208.root")'
```

```bash
~/venvs/v/bin/python ana/sim_psd.py output/tl208.root
```

```bash
~/venvs/v/bin/python ana/sim_psd.py output/bi214.root
```

```bash
root -l -b -q 'ana/sim_background.C("output/tl208.root=Tl208,output/bi214.root=Bi214")'
```

---

## ⓪ Geometry

### `geom/geo_check.mac` — read/build geometry

```
register the Ge detectors from the GDML's own map
overlap check = NOT {SKIP}, {NPOINTS} points per volume surface
load {GDML}; /run/initialize
```

```bash
remage --ignore-warnings -s GDML=$G -s SKIP=false -s NPOINTS=1000 -- geom/geo_check.mac 2>&1 | grep "Overlap is detected"
```

### `geom/geo_rt.C` — analysis on RT for cleanliness

```
[1] solids          the six material polycones bounding the wall
[2] clean surface   on-axis points = 2; duplicates, spikes, self-intersections = 0
[3] endcap flush    endcap radius (global max) vs barrel radius (mid-OFHC), within 0.05 mm
[4] shells          barrel radius vs the OFHC and SS outer bounds, within 0.05 mm
[5] wall thickness  r_wall - r_argon at 7 heights
[6] diff            extent, radii and seams against the mint geometry ../../l1000.gdml
[7] draw            whole tube, endcap, junction
```

---

## ① Simulation

### `sim/src_events.mac` — the source

```
seed {SEED}
confine uniformly by volume to reentrancetube + ofhc_cu + ss_316l
    reentrancetube = EFCu mother: endcap, lower wall AND the lid at the top
    ofhc_cu, ss_316l = the two shells placed inside it
ion {Z} {A} at rest, decaying that nuclide only
beamOn {NEV}
```

All three volumes, or only one section gets sampled. The wall is thin: remage needs ~550
tries per vertex, so the trial limit is raised to 100 000.

### `sim/sim_run.mac` — the run

```
register Ge from the GDML map, undergroundlar as a scintillator (the argon veto)
trees named by volume; gamma angular correlation on; overlap check off (geo_check did it)
load {GDML}; /run/initialize
store Ge and argon deposits with track ids, single precision; drop zero-energy Ge hits
store the gamma table
drop Ge, argon and gammas of every event with no Ge deposit (~99%); vertices are always kept
execute sim/src_events.mac
```

`Z`/`A` pick the chain, `NEV` the decays, `SEED` the random stream. Give every run its own
seed: the vertex is drawn first, so two chains with one seed decay at the same points and
their errors are no longer independent.

The event cut loses nothing `sim_background.C` uses (verified: identical output) and makes
files 12–17× smaller; the vertex table still counts every decay. 1 M decays take ~10 min at
`-t 8`. Output: `output/<run>.root` (Tl208 ~32 MB, Bi214 ~26 MB), directory `stp/` with one
tree per Ge detector (`V0101`…, 336), `undergroundlar`, `vtx`, `tracks`, `detector_origins`.

### `ana/src_vertices.C` — did the decays fill the tube

```
for each vertex: rt.materialAt(r, z) -> EFCu / OFHC / SS / LAr / outside
[1] read the vertices   [2] material split   [3] 16 slabs along z
verdict = nothing in argon, nothing outside, no empty slab
```

Expect `EFCu 11.8% / OFHC 43.3% / SS 44.9%`, `PASS`. A truly uniform fill would give EFCu
14.3% — see [Normalisation](#normalisation).

### `ana/sim_hits.C` — what a hit is (optional, for understanding)

```
[1] read every Ge step
[2] group by (decay, detector): steps -> responses -> responses above 5 keV -> decays
[3] energy share per particle
[4] the response > 100 keV with the most steps, step by step
[5] responses traced to the gamma that delivered most of their energy, with time and argon energy
[6] draw r-z map | x-y map | step vs summed spectrum | steps per response | depth | the [4] response
```

```bash
root -l -b -q 'ana/sim_hits.C("output/tl208.root")'
```

Output: `output/tl208_hits.png`. For Tl208: 71 706 steps → 6710 responses → 6128 decays;
electrons deposit 98.5% of the energy.

---

## ② Detector response — `ana/sim_psd.py`

Edgar's reboost `geds` group, operation for operation, with reboost's own functions.

```
per Ge detector, steps grouped per decay:
  active    = piecewise_linear_activeness(dist_to_surf, FCCD 1 mm, DLF 0.5)   dead layer
  energy    = sum(edep x active)
  drift     = 0° and 45° drift-time maps, blended by azimuth
  A_max     = maximum_current(...), smeared with the current resolution (seeded)
  AoE_class = (A_max/energy / mu - 1) / sigma(energy)
write evtid, det, energy_keV, aoe_class       (responses with energy > 0 only)
```

Parameters are Edgar's, written into the file; only the drift-time map is read from the
`Edgars_sim` mirror (reference detector `V00000A`, used for every detector). Reads the ROOT
file with uproot. A second argument changes the seed (default 1); one seed, one CSV.

Output: `output/<run>_psd.csv` (Tl208 6269 rows, Bi214 3126). Without it, `sim_background.C`
uses the raw deposit and stops at M1 + argon, and says so.

---

## ③ Background study — `ana/sim_background.C`

```
energy   active energy from <run>_psd.csv; a response with no row (all dead layer) is dropped
cuts     M1     exactly one detector above 5 keV in the decay
         argon  total undergroundlar deposit in the decay <= 20 keV   (where the 4 PE cut sits)
         PSD    AoE_class > -1.80
         ROI    |E - 2039| <= 55 keV
array    1000 kg Ge

input    list of run=isotope, wildcards allowed; each joined with output/<run>_psd.csv on
         (event, detector). files of one isotope merge into one run, events renumbered
memory   every decay is kept as counts; per-decay data only for decays with a Ge hit, so
         10^9-decay job arrays fit (1 M per chain: 0.9 GB, 10 s)

[1] configuration   [2] read one run   [3] combine the runs
per run
  [4] cut ladder       all hits: none -> M1 -> +argon -> +PSD
  [5] normalisation    per physical volume: V, M, A, decays in 10 yr, MC decays, MC density, weight
  [6] reach vs depth   decays, hits, hits per decay, ROI hits in 10 depth slabs
  [7] sizing           decays for 10% on the BI after cuts, and the time here; the dead zone
all runs summed
  [8] design table     ROI rate before / after cuts for 7 designs, against the goal
  [9] survival         AC, PSD, PSD|AC, Combined in the ROI, per chain x section
  [10] BI              per chain x section: BI +- radioassay +- MC statistics
  [11] draw            weighted spectrum (4 cut stages) | (L1, L2) seam scan | survival vs CDR and Edgar
```

Output: `output/Tl208_Bi214_background.png`, `output/Tl208_Bi214_survival.png`.

With the 1 M runs, as built, before cuts: BI `< 2.0e-06` (`0.20 ×` goal) — Tl208 `< 1.8e-06`
from 43 ROI hits (±15% MC; its EFCu activity is an upper limit), Bi214 `2.1e-07` from 2.
None of the 45 survive the cuts, so after cuts there is only a limit (`< 1e-07`).
**The as-built tube is under the goal before any cut.** Whether other materials would be
needs after-cut rates — see [Sizing](#sizing-a-run).

### Designs

A design is two seams, `L1` and `L2`, measured down from the top, and the material of the
three slabs they make. As built: steel to `L1 = 2.062 m`, Cu to `L2 = 4.067 m`, EFCu below
(the KS seams, read from the GDML). `[8]` compares:

```
KS as built: steel / Cu / EFCu        all EFCu     all Cu     all steel
KS, steel reaching twice as far down  KS, no steel (Cu down to L2)
same slabs, order reversed: EFCu / Cu / steel
```

The scan plot sweeps both seams (before cuts, as-built order, KS marked). All-Cu is `2.5 ×`
the goal before cuts: the cuts decide whether that design passes.

Reweighting changes a slab's activity, not the attenuation Geant4 already applied. That holds
while the wall thickness is unchanged; a thicker or thinner wall needs a new GDML and run.

### Normalisation

```
for each simulated decay i, at depth d_i, drawn in physical volume p_i (mother, OFHC, SS):
  m      = material the design puts at d_i
  real   = rho_m x A_m                        real decays per s per m^3 of that material
  MC     = n_p / V_p                          simulated decays per m^3, counted in p_i
  w_i    = real x 1 yr / MC / (1000 kg x 110 keV)     cts/(keV kg yr) for one simulated decay
BI = sum of w_i over ROI hits
```

**Why measure the MC density:** remage fills the mother (`reentrancetube`) ~20% more
sparsely than its two daughters; `[5]` prints the ratio (0.81 for Tl208, 0.80 for Bi214).
Assuming a uniform fill would put the BI ~23% low. Counting decays per volume makes the
weight right however remage samples. `V_p` comes from integrating the GDML polycones; an
independent point-in-solid integration agrees within 0.4%.

**Depth attribution:** the design assigns material by depth, so the EFCu lid at the top
counts as part of the steel slab. It changes nothing: no decay above ~3.1 m reaches a detector (`[6]`).

### Uncertainties

```
sigma_stat   sqrt(sum of w_i^2)               MC statistics
sigma_act    BI x dA/A                        radioassay
             printed separately; combined in quadrature they are survival_BI.py's error
no ROI hit   < 2.30 x mean weight of the decays there (90% CL), never zero; such a section
             stays out of the totals, which sum only what was measured
A is a limit the BI is a "<"
survival     binomial in ROI hits. PSD|AC = of the hits AC kept, how many PSD keeps too;
             Combined = AC x PSD|AC
```

Where sigma_act dominates, more simulation cannot help; only a better activity measurement can.

### Activities — `cfg::act`

| µBq/kg | Bi214       | Tl208   | source                              |
| ------ | ----------- | ------- | ----------------------------------- |
| steel  | 2500        | 1000    | Ralph's materialMix, no uncertainty |
| Cu     | 1           | 1       | Ralph's materialMix, no uncertainty |
| EFCu   | 0.19 ± 0.10 | < 0.077 | Edgar's `survival_BI.py` radioassay |

EFCu Tl208 is Edgar's upper limit (Ralph had 0.37). The PSD cut is Edgar's −1.80; his
`aoe_class_paras.yaml` separately lists −0.83. Survival references in the plot: CDR RE vessel
(EFCu) and Edgar's remage RE Cu.

### Sizing a run

`cfg::targetRel = 0.10`: a 10% BI **after cuts** needs 100 ROI hits that survive every cut.
`[7]` measures that rate once 10 hits survive; until then it estimates it as the before-cut
ROI rate × Edgar's Combined survival (Tl208 0.25%, Bi214 2.1%), and says so.

From the 1 M runs: **Tl208 9.3 × 10⁸ decays**, **Bi214 2.4 × 10⁹**. The Bi214 number rests on
2 ROI hits, so it is uncertain by ~70%: rerun ③ after part of the NERSC array and let `[7]`
correct it. Before cuts is already fine: 43 Tl208 ROI hits give 15%.

The dead zone: 71% of decays land above ~3.1 m and reach nothing. Confining the source below it
could buy up to 3–6×, but confining to a single volume stalled remage's sampler before (the
reason `src_events.mac` uses all three), so it is untested.

---

## Approximations vs Edgar's full chain

|            | here                            | Edgar                                            |
| ---------- | ------------------------------- | ------------------------------------------------ |
| argon veto | deposit ≤ 20 keV                | NPE ≥ 4 from the optical map                     |
| dead layer | distance to the nearest surface | distance to the n+ contact (differ only near p+) |
| PSD        | same reboost chain              | same                                             |

---

## NERSC — 10% after cuts

The sizing above, as job arrays of 10⁷ decays: **93 Tl208 tasks, 240 Bi214 tasks**, 32 threads
each on the `shared` QOS. ① runs in the same remage (the v0.26.0 container); ② and ③ run here,
where the Python and ROOT setup is tested. The task id is both the seed and the file name, so
every file is independent: more tasks later simply add statistics, and `scancel` drops extras.

**Here, once** — from `rt_sim/`, with your NERSC login and repo path filled in:

```bash
export NH=<user>@perlmutter.nersc.gov NDIR=<path>/LEGEND-background-simulation
```

```bash
ssh "$NH" "mkdir -p $NDIR/l1000_sim/rt_sim/output $NDIR/l1000_sim/rt_sim/sim"
```

```bash
rsync -av ../../KSendcap_l1kGeometry.gdml "$NH:$NDIR/"
```

```bash
rsync -av sim/ "$NH:$NDIR/l1000_sim/rt_sim/sim/"
```

**On Perlmutter** — from `l1000_sim/rt_sim/`:

```bash
shifterimg pull docker:legendexp/remage:v0.26.0
```

```bash
export G=../../KSendcap_l1kGeometry.gdml
```

```bash
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --array=1001-1093 --job-name=rt_tl208 --output=output/%x_%a.log --image=docker:legendexp/remage:v0.26.0 --wrap "shifter remage -q --ignore-warnings -t 32 -w -o output/tl208_\$SLURM_ARRAY_TASK_ID.root -s GDML=$G -s Z=81 -s A=208 -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -- sim/sim_run.mac"
```

```bash
sbatch --account=m2676 --constraint=cpu --qos=shared --time=04:00:00 --cpus-per-task=32 --array=2001-2240 --job-name=rt_bi214 --output=output/%x_%a.log --image=docker:legendexp/remage:v0.26.0 --wrap "shifter remage -q --ignore-warnings -t 32 -w -o output/bi214_\$SLURM_ARRAY_TASK_ID.root -s GDML=$G -s Z=83 -s A=214 -s NEV=10000000 -s SEED=\$SLURM_ARRAY_TASK_ID -- sim/sim_run.mac"
```

```bash
squeue --me
```

`$G` expands when you submit; `\$SLURM_ARRAY_TASK_ID` when the task runs. A task should take
about an hour (scaled from this laptop; the first log gives the real rate); 4 h is margin.
Output per task ~320 MB (Tl208) or ~260 MB (Bi214), ~90 GB in all, mostly the vertex table,
plus one `output/rt_*.log` per task.

**Here, after** — fetch the runs, then ② and ③ exactly as locally:

```bash
rsync -av --include='*_[0-9]*.root' --exclude='*' "$NH:$NDIR/l1000_sim/rt_sim/output/" output/
```

```bash
for f in output/tl208_*.root output/bi214_*.root; do ~/venvs/v/bin/python ana/sim_psd.py "$f"; done
```

```bash
root -l -b -q 'ana/sim_background.C("output/tl208_*.root=Tl208,output/bi214_*.root=Bi214")'
```

The wildcard expands inside the macro; each isotope's files merge into one run, events
renumbered, nothing counted twice. `[5]` confirms with `merged from 93 files` (Tl208).
At full statistics ② and ③ take ~3 h each, scaled from the 1 M runs.

---

## `geom/rt.h`

```
Outline                  (r, z) corners of one GDML polycone
  .radiusAt(z)           outermost radius at z
  .contains(r, z)        point inside?  (even-odd ray cast)
  .near / .deep(r, z, tol)  within tol of it / at least tol inside, tested in r AND z
  .areaAt(z)             cross-section area; rings handled (the argon at the lid)
  .volume()              stacked cross-sections
rtRead(file, solid)      one polycone out of the GDML text
RT = rtLoad(gdml)        wall, argon, OFHC and SS shell bounds, zBottom/zTop, seamOFHC/seamSS
  .sectionAt(z)          EFCu | OFHC | SS by height alone
  .materialAt(r, z)      EFCu | OFHC | SS | LAr | outside, exact (the lid is EFCu), 1 µm slack
                         in r and z for float32 vertices (0.2 µm in z is 40 µm in r on the lid cone)
  .physAt(r, z)          physical volume the source drew from: 0 mother, 1 OFHC, 2 SS
  .physVolume(p)         its volume
  .wallVolume(z0, z1)    wall volume between two heights
rtResolve(gdml)          which GDML: argument -> $RT_GDML -> search upward
rtIsGermanium(name)      is this tree a Ge detector ("V" + 4 digits)?
rtScan(tree, cols, f)    walk a tree row by row, handing the named columns to f
rtOut(name)              "output/<name>", creating output/
rtVerdict(tag, ok)       PASS / FAIL line and batch exit code
```
