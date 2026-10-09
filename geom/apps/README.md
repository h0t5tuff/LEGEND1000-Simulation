# apps

Double-click apps for the two autopeels, each with an icon of its own detector array: `/Applications/l200.app`
peels `LEGEND200/`'s geometry, `/Applications/l1000.app` peels `geom/`'s. Standalone, an app carries its geometry
and a native viewer, so it runs on any Mac with nothing installed; as a launcher, it runs the repo's autopeel.

| file | |
| ---- | --- |
| `make_app.py` | builds one app: its icon, the standalone viewer + geometry (or the launcher), the bundle; installs it |
| `autopeel.swift` | the standalone viewer: the autopeel in Swift/SceneKit, reading `geometry.json` + `geometry.bin` from its bundle |
| `output/<name>/render.png` | the HPGe detectors, string hardware, ultem and PEN of the autopeel, rendered offscreen (1400 px, transparent) |
| `output/<name>/art.png` | the render above the label (`L200`, `L1000`) on a transparent 1024 px canvas |
| `output/<name>/AppIcon.icon` | the Icon Composer package: `icon.json` (the background gradient, the art as one plain layer) and `Assets/art.png` |
| `output/<name>/<name>.app` | the standalone app, as installed in `/Applications` |
| `output/<name>/<name>.zip` | the same, zipped, to give to someone |

```bash
~/venvs/v/bin/python geom/apps/make_app.py l200
~/venvs/v/bin/python geom/apps/make_app.py l1000
~/venvs/v/bin/python geom/apps/make_app.py l200 --launcher
/Applications/l200.app/Contents/MacOS/l200 --snapshot out.png --peel 12
```

- **Standalone (the default):** `Contents/MacOS/<name>` is `autopeel.swift` compiled for Apple silicon and Intel
  (one universal binary, macOS 12 or later, only system frameworks); `Contents/Resources/geometry.json` holds the
  layers (the autopeel's `LAYERS`, in its order and with its labels), each split by colour and hierarchy depth into
  groups of parts, and `geometry.bin` each part's mesh once plus all its placements, so the app builds the full
  scene when it opens. l200 on 8 Oct 2026: 7 s to build, 3.8 MB of geometry, a 6.4 MB app, a 3.1 MB zip; it opens
  in 0.6 s with 1.59 M triangles in 20 layers and uses ~240 MB (the Python autopeel ~540 MB); l1000: 31 s (its
  autopeel re-meshed first), 8.4 MB of geometry, an 11 MB app, a 3.4 MB zip, 2.80 M triangles in 19 layers, ~400 MB. It peels like the
  autopeel (3 s per layer, 0.5 s fades, the same keys and text); drag rotates around the vertical axis, scroll or
  pinch zooms (SceneKit's camera control). A daughter sharing a face with its mother is drawn a hair towards the eye, as the autopeel's polygon offset.
  The geometry is what the autopeel had when the app was built: rebuild after the GDML changes.
- **`--snapshot out.png [--peel N] [--size WxH]`:** renders one frame offscreen, N layers peeled and re-framed, and
  quits; the check that the scene is right without opening the window.
- **`--launcher`:** `Contents/MacOS/<name>` is a bash script that runs the repo's autopeel with `~/venvs/v`, its
  printout in `$TMPDIR/<name>-autopeel.log` and an alert if it fails; it holds this repo's absolute path.
- **The icon:** macOS 27 shows a plain `.icns` icon (and a classic Xcode icon set) shrunk inside a dark frame; an
  `AppIcon.icon` compiled into `Assets.car` by `actool`, named by `CFBundleIconName`, is drawn as a native icon.
  Its gradient is set for the dark appearance too, else dark icon mode swaps in a black background. Change the
  camera (`view`), colours (`fill`) or label in `APPS`, the layers drawn in `SHOW`, and rerun; hand edits to
  `AppIcon.icon` in Icon Composer are overwritten by the next run.
- **Giving it to someone:** send `<name>.zip`. The app is signed ad hoc, which runs anywhere but is no identity
  Apple vouches for: copied by USB drive or `scp` it opens with a double-click; downloaded, mailed or AirDropped,
  macOS blocks the first open, and the recipient allows it once in System Settings → Privacy & Security → Open
  Anyway. No warning at all needs a Developer ID signature and Apple's notarization (Apple Developer Program,
  99 USD a year); the App Store needs the same membership, the app sandbox, and Apple's review.
- **Who may have it:** l200's geometry comes from the private legend-metadata. Share it within LEGEND; anything
  public (the App Store, a website) needs the collaboration's approval, or a geometry built with `--public-geom`.
- A running app is safe to rebuild: the copy in `/Applications` is replaced, and the new one starts next time.
  If the Dock keeps an old icon, it updates on its next refresh (or the next login).
