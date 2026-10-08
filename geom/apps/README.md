# apps

Double-click apps for the two autopeels: `/Applications/l200.app` runs `LEGEND200/l200-autopeel.py`,
`/Applications/l1000.app` runs `geom/l1000-autopeel.py`. Each has an icon of its own detector array.

| file | |
| ---- | --- |
| `make_app.py` | builds one app: renders its icon, packs it for Icon Composer, writes the bundle in `/Applications` |
| `output/<name>/render.png` | the HPGe detectors, string hardware, ultem and PEN of the autopeel, rendered offscreen (1400 px, transparent) |
| `output/<name>/art.png` | the render above the label (`L200`, `L1000`) on a transparent 1024 px canvas |
| `output/<name>/AppIcon.icon` | the Icon Composer package: `icon.json` (the background gradient, the art as one plain layer) and `Assets/art.png` |

```bash
~/venvs/v/bin/python geom/apps/make_app.py l200
~/venvs/v/bin/python geom/apps/make_app.py l1000
```

- **What it builds:** `/Applications/<name>.app/Contents/` holds `Info.plist`, `MacOS/<name>` (the bash launcher:
  it runs the autopeel with `~/venvs/v`, writes its printout to `$TMPDIR/<name>-autopeel.log` and shows an alert
  if it fails) and `Resources/Assets.car` + `AppIcon.icns`, compiled from `AppIcon.icon` by Xcode's `actool`.
  3 s per app (8 Oct 2026), plus ~20 s the first time if the autopeel has never meshed its geometry.
- **Why Icon Composer:** macOS 27 shows a plain `.icns` icon (and a classic Xcode icon set) shrunk inside a
  dark frame. An `AppIcon.icon` compiled into `Assets.car`, named by `CFBundleIconName`, is drawn as a native
  icon. Its gradient is set for the dark appearance too, else dark icon mode swaps in a black background.
- **Changing an icon:** edit `APPS` (the camera `view`, the colours `fill`, the label) or `SHOW` (the layers
  drawn) and rerun. Hand edits to `output/<name>/AppIcon.icon` in Icon Composer (Xcode → Open Developer Tool)
  are overwritten by the next run; to keep them, compile that package with the `actool` call in `build_app`.
- **The launchers hold this repo's absolute path:** after moving the repo, rerun `make_app.py`. A running
  app is safe to rebuild: the launcher is replaced, not edited in place.
- If the Dock keeps an old icon, it updates on its next refresh (or the next login).
