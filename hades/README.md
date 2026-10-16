# hades

HPGe characterization at the HADES test stand: calibration sources fired at one detector, to measure
its dead layer (FCCD), active volume and resolution. The stand comes from
[legend-pygeom-hades](https://legend-pygeom-hades.readthedocs.io); run from this folder.

| source | measures | observable |
| ------ | -------- | ---------- |
| Am-241 (`am_HS1`, collimated) | dead layer | 59.5 keV peak per decay |
| Ba-133 (`ba_HS4`) | dead layer | R = (79.6+81.0) / 356.0 keV peaks, no activity needed |
| Co-60 (`co_HS5`) | active volume | 1332.5 keV peak per decay |
| Th-228 (`th_HS2`) | resolution | the chain's peaks from 239 to 2615 keV |

`sim/<src>.yaml` builds each stand, `sim/run.mac` decays the source (Am-241: its 59.5 keV line aimed
into a 10° cone) and records the steps in the crystal, and `ana/spectrum.C` applies the dead layer
and resolution, counts the peaks and scans the FCCD. Everything generated goes to `output/`.

```bash
~/venvs/v/bin/legend-pygeom-hades --public-geom --config sim/co60.yaml output/co60.gdml
remage -q --ignore-warnings -t 8 -w -o output/co60.root -s GDML=output/co60.gdml -s PARTICLE=ion -s Z=27 -s A=60 -s E=0 -s THETA=180 -s NEV=2000000 -s SEED=3 -- sim/run.mac
root -l -b -q 'ana/spectrum.C("output/co60.root")'
```

Same for `ba133` (`Z=56 A=133 SEED=4`), `th228` (`Z=90 A=228 NEV=5000000 SEED=2`) and `am241`
(`PARTICLE=gamma Z=95 A=241 E=59.5409 THETA=10 NEV=4000000 SEED=1`): one seed per run. A measured observable reads back as
an FCCD: `ana/spectrum.C("output/ba133.root", 1.0, ‹R›, ‹error›)`.

Results (2-5 × 10⁶ primaries, public dummy detector, FCCD 1 mm; remage 0.26, 6 Oct 2026):

| | Am-241 | Ba-133 | Co-60 | Th-228 |
| --- | --- | --- | --- | --- |
| energy lost to the dead layer | 54% | 20% | 7.7% | 8.5% |
| observable | 7.69 × 10⁻⁶ /decay | R = 0.799 | 9.82 × 10⁻³ /decay | 1.48 × 10⁻³ /decay |
| a 1% error on it moves the FCCD by | 0.010 mm | 0.026 mm | 0.051 mm | 0.050 mm |

Limits: an effective FCCD (distance to the nearest surface, linear transition), a dummy detector, a
placeholder resolution, MC errors only. Future work: ‹field map / A/E› · ‹pulse shapes› · ‹real
detector› · ‹n+ contact› · ‹transition layer› · ‹resolution calibration› · ‹activity and live time› ·
‹surface scans› · ‹DAQ effects› · ‹systematics›.
