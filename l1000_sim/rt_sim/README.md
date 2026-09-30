# RT backgrounds — KS-endcap L1000

**Question.** How much background does the re-entrant tube (RT) put into the Ge detectors,
and does the choice and placement of its materials matter?
**Target.** `1e-5 cts/(keV·kg·yr)` in the ROI, `2039 ± 55 keV`.

**Idea.** Simulate each decay chain once, uniformly through the whole tube wall, then give
every decay a weight for the material a design puts at its depth. Any design — which
material, where the seams sit — becomes arithmetic on one run, not a new simulation.

---

## Pipeline

```
⓪ design     geo_check.mac   does the KS GDML build, does anything overlap
             geo_rt.C        is the tube one clean surface               -> output/geo_rt.png
① remage     sim_run.mac     + src_events.mac: decays in the tube wall   -> output/<run>.root
             src_vertices.C  did the decays fill the whole tube
             sim_hits.C      what a Ge hit is made of                    -> output/<run>_hits.png
② response   sim_psd.py      dead layer + A/E, via reboost               -> output/<run>_psd.csv
③ cuts       sim_background.C  cuts, normalisation, designs, BI          -> output/<Isos>_background.png
                                                                            output/<Isos>_survival.png
```

One run per chain: **Tl208** (Th232) and **Bi214** (U238). The background is their sum.
Every generated file lands in `output/`. `rt_geom.h` is the shared geometry reader:
no macro hard-codes a radius or a seam, every number comes from the GDML.

## Setup

From `l1000_sim/rt_sim/`:

```bash
export PATH=/Users/tensor/Documents/REMAGE/install-remage-v0.25.0/bin:$PATH
```

```bash
export G=../../KSendcap_l1kGeometry.gdml
```

remage v0.25.0 (Geant4 11.4), ROOT 6.40, Python in `~/venvs/v` (reboost 1.3.1, uproot).
The `.C` macros find the GDML themselves (argument → `$RT_GDML` → search upward).

## Run everything

```bash
remage -q --ignore-warnings -s GDML=$G -s SKIP=true -s NPOINTS=1000 -- geo_check.mac
```
```bash
root -l -b -q geo_rt.C
```
```bash
remage -q --ignore-warnings -t 8 -w -o output/tl208.root -s GDML=$G -s Z=81 -s A=208 -s NEV=100000 -s SEED=1 -- sim_run.mac
```
```bash
remage -q --ignore-warnings -t 8 -w -o output/bi214.root -s GDML=$G -s Z=83 -s A=214 -s NEV=100000 -s SEED=1 -- sim_run.mac
```
```bash
root -l -b -q 'src_vertices.C("output/tl208.root")'
```
```bash
~/venvs/v/bin/python sim_psd.py output/tl208.root
```
```bash
~/venvs/v/bin/python sim_psd.py output/bi214.root
```
```bash
root -l -b -q 'sim_background.C("output/tl208.root=Tl208,output/bi214.root=Bi214")'
```

Each step is explained below.

---

## ⓪ Geometry

### `geo_check.mac` — does it build

```
register the Ge detectors from the GDML's own map
overlap check = NOT {SKIP}, {NPOINTS} points per surface
load {GDML}; /run/initialize
```

`SKIP=true` builds only (~6 s). The overlap scan (~3 min):

```bash
remage --ignore-warnings -s GDML=$G -s SKIP=false -s NPOINTS=1000 -- geo_check.mac 2>&1 | grep "Overlap is detected"
```

### `geo_rt.C` — is the tube clean

```
[a] the six solids bounding the wall
[b] on-axis points = 2   [c] duplicate points = 0   [d] spikes = 0   [e] self-intersections = 0
[f] endcap radius vs barrel radius, within 0.05 mm
[g] barrel vs the OFHC and SS shell outer bounds
[h] wall thickness at 7 heights: 1.5 mm EFCu, 6 mm OFHC and SS
[i] every number diffed against the mint geometry ../../l1000.gdml
verdict = clean surface AND endcap flush AND shells agree
```

Output: `output/geo_rt.png` (whole tube, endcap, junction) and `RESULT [geo_rt]`.
It reports **FAIL** today, for one real reason: the KS endcap is `1.0 mm` wider than the
barrel at `z = −607 mm`. Everything else passes.

---

## ① Simulation

### `src_events.mac` — the source

```
confine uniformly by volume to reentrancetube + ofhc_cu + ss_316l
    reentrancetube = EFCu mother: endcap, lower wall AND the lid at the top
    ofhc_cu, ss_316l = the two shells placed inside it
ion {Z} {A} at rest, decaying that nuclide only
beamOn {NEV}
```

All three volumes, or only one section gets sampled. The wall is thin (~560 tries per
vertex), so the trial limit is raised to 100 000.

### `sim_run.mac` — the run

```
register Ge from the GDML map, undergroundlar as a scintillator (the argon veto)
store Ge + argon deposits with track ids, and the gamma table
gamma angular correlation on; overlap check off (geo_check.mac did it)
execute src_events.mac
```

`Z`/`A` pick the chain, `SEED` the random stream, `NEV` the decays.
Output: `output/<run>.root`, trees `stp/V*` (Ge), `stp/undergroundlar`, `stp/vtx`, tracks.

### `src_vertices.C` — did the decays fill the tube

```
for each vertex: rt.materialAt(r, z) -> EFCu / OFHC / SS / LAr / outside
[a] material split   [b] 16 slabs along z
verdict = nothing in argon, nothing outside, no empty slab
```

Expect `EFCu 12% / OFHC 43% / SS 45%`, then `PASS`. EFCu would be 14% if the fill were
truly uniform by volume — see [Normalisation](#normalisation).

### `sim_hits.C` — what a hit is (optional, for understanding)

```
read every Ge step; group by (decay, detector) into one response
[a] steps -> responses -> decays      [b] energy share per particle
[c] one response, step by step        [d] each response traced back to its gamma
```

```bash
root -l -b -q 'sim_hits.C("output/tl208.root")'
```

Output: `output/tl208_hits.png`.

---

## ② Detector response — `sim_psd.py`

Edgar's reboost `geds` group, operation for operation, with reboost's own functions.

```
per Ge detector, steps grouped per decay:
  active    = piecewise_linear_activeness(dist_to_surf, FCCD 1 mm, DLF 0.5)   dead layer
  energy    = sum(edep x active)
  drift     = 0° and 45° drift-time maps, blended by azimuth
  A_max     = maximum_current(...), smeared with the current resolution (seeded)
  AoE_class = (A_max/energy / mu - 1) / sigma(energy)
write evtid, det, energy_keV, aoe_class
```

Parameters are Edgar's, written into the file; only the drift-time map is read from the
`Edgars_sim` mirror (reference detector `V00000A`). Reads the ROOT file with uproot.
A second argument changes the seed (default 1); the same seed gives the same CSV.

Output: `output/<run>_psd.csv`. Without it, `sim_background.C` stops at M1 + argon and says so.

---

## ③ Background study — `sim_background.C`

```
cuts:  M1        only one detector above 5 keV
       argon     argon deposit <= 20 keV            (where the 4 PE cut sits)
       PSD       AoE_class > -1.80                   (from the _psd.csv)
       ROI       2039 +- 55 keV, on the active energy
array: 1000 kg Ge, 10 yr

input: list of run=isotope, joined with output/<run>_psd.csv on (event, detector)
       several files of one isotope are merged into one run (job arrays)

per run
  [a] cut ladder       hits: none -> M1 -> +argon -> +PSD
  [b] normalisation    per physical volume: V, M, A, N_real, MC decays, MC density, weight
  [c] reach vs depth   hits per decay in 10 depth slabs
  [d] sizing           decays needed for 10% / 5%; the dead zone and the speed-up from skipping it
all runs summed
  [e] design table     ROI rate before / after cuts for 7 designs, against the goal
  [f] survival         AC, PSD, PSD|AC, Combined in the ROI, per chain x section
  [g] BI               per chain x section: BI +- radioassay +- MC statistics
plots: weighted spectrum (4 cut stages) | (L1, L2) seam scan | survival vs CDR and Edgar
```

Output: `output/Tl208_Bi214_background.png`, `output/Tl208_Bi214_survival.png`.

### Designs

A design is two seams, `L1` and `L2`, measured down from the top, and the material of the
three slabs they make. As built: steel to `L1 = 2.062 m`, Cu to `L2 = 4.067 m`, EFCu below
(the KS seams, read from the GDML). `[e]` compares all-EFCu, all-Cu, all-steel, moved seams
and the reversed order; the scan plot sweeps both seams.

Reweighting changes a slab's activity, not the attenuation Geant4 already applied. That holds
while the wall thickness is unchanged; a thicker or thinner wall needs a new GDML and run.

### Normalisation

```
for each simulated decay i, at depth d_i, drawn in physical volume p_i:
  m        = material the design puts at d_i
  real     = rho_m x A_m x T               real decays per m^3 of that material in T
  MC       = n_p / V_p                     simulated decays per m^3, counted in p_i
  w_i      = real / MC / (M_Ge x 110 keV)  cts/(keV kg yr) for one simulated decay
BI = sum of w_i over ROI hits
```

**Why measure the MC density:** remage fills the mother (`reentrancetube`) about 19% more
sparsely than its two daughters (`[b]` prints the ratio, 0.814). Assuming a uniform fill put
the BI ~23% low. Counting decays per volume makes the weight correct however remage samples.
`V_p` comes from integrating the GDML polycones; an independent point-in-solid integration
agrees within 0.4%.

**Depth attribution:** the design assigns material by depth, so the EFCu lid at the top
counts as part of the steel slab. It changes nothing: no decay above ~3 m reaches a detector (`[c]`).

### Uncertainties

| | |
| --- | --- |
| BI | `sqrt(σ_act² + σ_stat²)`, as `survival_BI.py` |
| σ_act | `BI × ΔA/A`, the radioassay |
| σ_stat | `sqrt(Σ w_i²)`, the MC |
| no ROI hit | `< 2.30 × mean weight` (90% CL), never zero |
| activity is a limit | the BI is a `<` |
| survival | binomial in ROI hits; `PSD|AC` = of the hits AC kept, how many PSD keeps; `Combined = AC × PSD|AC` |

Where σ_act dominates, more simulation cannot help; only a better activity measurement can.

### Activities — `cfg::act`

| µBq/kg | Bi214 | Tl208 | source |
| --- | --- | --- | --- |
| steel | 2500 | 1000 | Ralph's materialMix, no uncertainty |
| Cu | 1 | 1 | Ralph's materialMix, no uncertainty |
| EFCu | 0.19 ± 0.10 | < 0.077 | Edgar's radioassay |

EFCu Tl208 is Edgar's upper limit (Ralph had 0.37). The PSD cut is Edgar's −1.80; his
`aoe_class_paras.yaml` separately lists −0.83.

### Sizing a run

Relative error `p` needs `1/p²` ROI hits, so `N = 1/(p² × ROI hits per decay)`. `[d]` prints
this from the run just made, and the dead zone: the depth above which no decay reaches a
detector (~3 m, ~70% of all decays). Confining the source below it would buy ~3–6× statistics.

---

## Approximations vs Edgar's full chain

| | here | Edgar |
| --- | --- | --- |
| argon veto | deposit ≤ 20 keV | NPE ≥ 4 from the optical map |
| dead layer | distance to the nearest surface | distance to the n+ contact (differ only near p+) |
| PSD | same reboost chain | same |

---

## NERSC

Write `.root`: `sim_psd.py` and `sim_background.C` both read it. One job per array task,
each with its own `SEED`:

```bash
remage -q --ignore-warnings -t "$SLURM_CPUS_PER_TASK" -w -o "output/tl208_$SLURM_ARRAY_TASK_ID.root" -s GDML=$G -s Z=81 -s A=208 -s NEV=1000000 -s SEED="$SLURM_ARRAY_TASK_ID" -- sim_run.mac
```

Same for Bi214 (`Z=83 A=214`, `output/bi214_...`). Then:

```bash
for f in output/tl208_*.root output/bi214_*.root; do ~/venvs/v/bin/python sim_psd.py "$f"; done
```
```bash
root -l -b -q 'sim_background.C("output/tl208_*.root=Tl208,output/bi214_*.root=Bi214")'
```

The wildcard expands inside the macro; the files of each isotope merge into one run, so
nothing is counted twice.

---

## `rt_geom.h`

```
Outline                  (r, z) corners of one GDML polycone
  .radiusAt(z)           outermost radius at z
  .contains(r, z, tol)   point inside?  (ray cast)
  .areaAt(z)             cross-section area; rings handled (the argon at the lid)
  .volume()              stacked cross-sections
RT = rtLoad(gdml)        wall, argon, OFHC and SS shell bounds, zBottom/zTop, seamOFHC/seamSS
  .materialAt(r, z)      EFCu | OFHC | SS | LAr | outside
  .physAt(r, z)          physical volume the source drew from: 0 mother, 1 OFHC, 2 SS
  .physVolume(p)         its volume
  .wallVolume(z0, z1)    wall volume between two heights
rtOut(name)              "output/<name>", creating output/
rtVerdict(tag, ok)       PASS / FAIL line and batch exit code
```
