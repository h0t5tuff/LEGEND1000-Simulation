|     | step                                | result                                                 |
| --- | ----------------------------------- | ------------------------------------------------------ |
| 1   | geometry builds / overlap scan      | **PASS** with one known overlap, outside the RT        |
| 2   | `rt_verify.C`                       | **FAIL** a real 1 mm endcap step                       |
| 3   | `src_vertices.C`                    | **PASS**                                               |
| 4   | 1M-event production + `sim_reach.C` | **PASS** — reach is 0.4366 %, and 98.9 % of it is EFCu |

## 1. Geometry build and overlap scan

Two `Overlap is detected` lines, which are **one physical overlap reported from both
sides**:

```
skirt:0 (G4Tubs) with outercryostat:0 (G4GenericPolycone)
  overlap at local point (2020.76, 2674.4, -2679.44) by 5.04412 mm  (max of 5 cases)
outercryostat:0 (G4GenericPolycone) with skirt:0 (G4Tubs)
  overlap at local point (-875.786, 3238.58, 1260.31) by 2.91165 mm  (max of 5 cases)
```

## 2. Reentrance-tube surface

RT spans z −1259.0 … 4987.0 mm (length 6246.0), seams at 920.0 (EFCu\|OFHC) and
2925.0 (OFHC\|SS).

**the one failure.** The endcap head peaks at r = 999.9998 mm against a barrel
radius of 999.0000 mm, a **+0.9998 mm step** at z ≈ −607. `rt_verify.png` shows it
between z ≈ −560 and −690, then returns to 999.

**wall thickness** is sane throughout — 1.5 mm through the EFCu endcap and lower
wall, thickening to 3.0 mm at the seam approach, then 6.0 mm across both OFHC and SS.

**KS vs mint l1000** — the KS file is the mint geometry displaced by a uniform
**+68 mm**: both z ends, both seams, and the barrel radius all move +68.00, and the
length is preserved exactly (6246.00 both). The endcap rmax moves +68.9998, which is
the same +68 plus the 1 mm step.

## 3. Source vertices along the tube

5000 events in wall. **PASS.**

| section | vertices | share  |
| ------- | -------- | ------ |
| EFCu    | 601      | 12.0 % |
| OFHC    | 2170     | 43.4 % |
| SS      | 2229     | 44.6 % |

Vertices span z −1228.0 … 4987.0 mm, r 50.3 … 999.9 mm.

## 4. Production run — 1M events

**Exit 0 in 8 min 44 s** (3829 s CPU, 732 % — at this event count the 8 threads
saturate). Output `rt_ge.root`,

**4366 of 1,000,000 events reached an HPGe — 0.4366 ± 0.0066 %.**

### [a] By section

| section | fired   | reached | reach [%]                    |
| ------- | ------- | ------- | ---------------------------- |
| EFCu    | 118 316 | 4316    | **3.6479 ± 0.0545**          |
| OFHC    | 432 872 | 50      | 0.0116 ± 0.0016              |
| SS      | 448 812 | 0       | **0** (< 0.00067 at 95 % CL) |

**The EFCu endcap dominates completely: 98.9 % of every event
that reaches a detector starts there**, from only 12 % of the wall volume. Its
per-decay reach is **315× the OFHC section**.

The SS section produced zero hits from
448 812 decays despite being the single largest volume. SS is at least ~5000× quieter than EFCu.

consequence: RT background from this geometry is an **EFCu radiopurity
problem**

### [b] Reach vs height

| z [mm]      | section  | fired        | reached | reach [%]  |
| ----------- | -------- | ------------ | ------- | ---------- |
| −1129       | EFCu     | 4767         | 315     | 6.6079     |
| −869        | EFCu     | 10 132       | 660     | 6.5140     |
| −608        | EFCu     | 11 359       | 864     | **7.6063** |
| −348        | EFCu     | 11 375       | 838     | 7.3670     |
| −88         | EFCu     | 11 475       | 816     | 7.1111     |
| 172         | EFCu     | 11 327       | 531     | 4.6879     |
| 433         | EFCu     | 11 435       | 210     | 1.8365     |
| 693         | EFCu     | 11 419       | 61      | 0.5342     |
| 953         | OFHC     | 39 668       | 41      | 0.1034     |
| 1213        | OFHC     | 56 371       | 20      | 0.0355     |
| 1474        | OFHC     | 56 210       | 7       | 0.0125     |
| 1734        | OFHC     | 56 278       | 3       | 0.0053     |
| 1994 … 4857 | OFHC, SS | ~56 000 each | 0       | 0          |
