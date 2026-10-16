# app

`/Applications/LEGEND.app` opens with three viewers to choose from:

- **l200**: the LEGEND-200 geometry, peeled layer by layer as you rotate it
- **l1000**: the same for LEGEND-1000
- **rt**: the re-entrance tube designs: two sliders move the seams, and the 3D tube takes the colour of that
  design's background index

```bash
~/venvs/v/bin/python app/make_app.py
```

builds it, installs it and keeps `output/LEGEND.zip` to share (~30 s). Rebuild after a GDML or `output/designs.json`
changes.

| file             |                                         |
| ---------------- | --------------------------------------- |
| `make_app.py`    | the build                               |
| `legend.swift`   | the choices, and what the viewers share |
| `autopeel.swift` | l200 and l1000                          |
| `rt.swift`       | rt                                      |

Share the zip within LEGEND only (l200's geometry is private). A downloaded copy opens after one "Open Anyway"
in System Settings → Privacy & Security.
