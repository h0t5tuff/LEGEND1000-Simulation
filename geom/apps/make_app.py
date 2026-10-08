#!/usr/bin/env python
"""Build a double-click app in /Applications that runs one of the repo's autopeels, with an icon of its geometry.

    ~/venvs/v/bin/python geom/apps/make_app.py l200   (or l1000; see README.md)

1. render the detector array offscreen from the autopeel's geometry   -> output/<name>/render.png
2. put the label under it, on a transparent 1024 px canvas            -> output/<name>/art.png
3. wrap it in an Icon Composer package with the background gradient   -> output/<name>/AppIcon.icon
4. write /Applications/<name>.app: Info.plist, the launcher script, and the icon compiled by Xcode's actool
"""

import importlib.util
import json
import os
import plistlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

REPO = Path(__file__).resolve().parents[2]  # geom/apps/make_app.py -> the repository
OUT = Path(__file__).resolve().parent / "output"
APPS = {  # autopeel, title, label, camera (theta, phi from +z, deg; zoom), background gradient (top, bottom)
    "l200": {"peel": "LEGEND200/l200-autopeel.py", "title": "LEGEND-200", "label": "L200", "view": (62, 25, 1.0), "fill": ("#2b4f86", "#0d1830")},
    "l1000": {"peel": "geom/l1000-autopeel.py", "title": "LEGEND-1000", "label": "L1000", "view": (40, 30, 1.15), "fill": ("#1d6b6f", "#08262a")},
}
SHOW = ("string hardware", "ultem", "PEN")  # autopeel layers drawn with the HPGe detectors
FONT = "/System/Library/Fonts/SFNS.ttf"
LSREGISTER = "/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister"


def render(peel: Path, png: Path, theta: float, phi: float, zoom: float) -> None:
    """The HPGe core plus the SHOW layers of an autopeel, offscreen, on a transparent background."""
    from vtkmodules.vtkIOImage import vtkPNGWriter
    from vtkmodules.vtkRenderingCore import vtkRenderer, vtkRenderWindow, vtkWindowToImageFilter

    spec = importlib.util.spec_from_file_location("peel", peel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    layers, _ = m.assemble(m.load_geometry(), m.read_colours())

    ren = vtkRenderer()
    ren.SetBackgroundAlpha(0.0)
    for j, lay in enumerate(layers):
        if j == len(layers) - 1 or any(k in lay["label"] for k in SHOW):
            for actor, _ in lay["actors"]:
                p = actor.GetProperty()
                p.SetAmbient(0.25)
                p.SetSpecular(0.35)
                p.SetSpecularPower(30)
                ren.AddActor(actor)
    th, ph = np.radians(theta), np.radians(phi)
    cam = ren.GetActiveCamera()
    cam.SetFocalPoint(0, 0, 0)
    cam.SetPosition(np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th))
    cam.SetViewUp(0, 0, 1)
    ren.ResetCamera()
    cam.Zoom(zoom)
    ren.ResetCameraClippingRange()

    win = vtkRenderWindow()
    win.SetOffScreenRendering(1)
    win.SetAlphaBitPlanes(1)
    win.SetMultiSamples(8)
    win.SetSize(1400, 1400)
    win.AddRenderer(ren)
    win.Render()
    grab = vtkWindowToImageFilter()
    grab.SetInput(win)
    grab.SetInputBufferTypeToRGBA()
    grab.Update()
    w = vtkPNGWriter()
    w.SetFileName(str(png))
    w.SetInputConnection(grab.GetOutputPort())
    w.Write()


def art(render_png: Path, png: Path, label: str) -> None:
    """The render in the upper part of a transparent 1024 px canvas, the label in bold white below it."""
    S = 1024
    canvas = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    im = Image.open(render_png).convert("RGBA")
    im = im.crop(im.getchannel("A").getbbox())
    k = min(746 / im.width, 677 / im.height)  # fits a 746 x 677 box centred at y = 416
    im = im.resize((round(im.width * k), round(im.height * k)), Image.LANCZOS)
    canvas.alpha_composite(im, ((S - im.width) // 2, round(416 - im.height / 2)))
    font = ImageFont.truetype(FONT, 149)
    font.set_variation_by_name("Bold")
    d = ImageDraw.Draw(canvas)
    x0, y0, x1, y1 = d.textbbox((0, 0), label, font=font)
    x, y = (S - (x1 - x0)) // 2 - x0, 860 - (y1 - y0) // 2 - y0
    d.text((x, y + 4), label, font=font, fill=(0, 0, 0, 110))  # soft drop
    d.text((x, y), label, font=font, fill=(255, 255, 255, 255))
    canvas.save(png)


def icon_package(art_png: Path, pkg: Path, top: str, bottom: str) -> None:
    """An Icon Composer package: the gradient as fill (in the dark appearance too) and the art as one plain layer."""
    p3 = lambda h: "display-p3:" + ",".join(f"{int(h[i:i + 2], 16) / 255:.5f}" for i in (1, 3, 5)) + ",1.00000"
    shutil.rmtree(pkg, ignore_errors=True)
    (pkg / "Assets").mkdir(parents=True)
    shutil.copy(art_png, pkg / "Assets" / "art.png")
    gradient = {"linear-gradient": [p3(top), p3(bottom)]}
    icon = {
        "fill-specializations": [{"value": gradient}, {"appearance": "dark", "value": gradient}],
        "groups": [{
            "layers": [{"image-name": "art.png", "name": "art", "glass": False,
                        "position": {"scale": 1, "translation-in-points": [0, 0]}}],
            "shadow": {"kind": "neutral", "opacity": 0.35},
            "translucency": {"enabled": False, "value": 0.5},
            "specular": False,
        }],
        "supported-platforms": {"circles": ["watchOS"], "squares": "shared"},
    }
    (pkg / "icon.json").write_text(json.dumps(icon, indent=2) + "\n")


LAUNCHER = """#!/bin/bash
# {name}: double-click to peel the {title} geometry, i.e. run {peel}.
# Written by geom/apps/make_app.py. Its printout goes to $TMPDIR/{stem}.log; if it fails, the last lines show in an alert.
DIR="{dir}"
LOG="${{TMPDIR:-/tmp}}/{stem}.log"
export PYTHONUNBUFFERED=1  # the log fills while the window is open
alert() {{ osascript -e "display alert \\"{name}\\" message \\"$(printf '%s' "$1" | tr -d '"\\\\')\\"" >/dev/null; }}
cd "$DIR" 2>/dev/null || {{ alert "No folder $DIR"; exit 1; }}
"{python}" {script} >"$LOG" 2>&1 || alert "The autopeel stopped with an error:
$(tail -n 6 "$LOG")"
"""


def build_app(name: str, cfg: dict, pkg: Path) -> Path:
    app = Path("/Applications") / f"{name}.app"
    contents, peel = app / "Contents", REPO / cfg["peel"]
    (contents / "MacOS").mkdir(parents=True, exist_ok=True)
    (contents / "Resources").mkdir(exist_ok=True)
    info = {
        "CFBundleExecutable": name, "CFBundleIdentifier": f"local.legend.{name}", "CFBundleName": name,
        "CFBundleDisplayName": name, "CFBundlePackageType": "APPL", "CFBundleShortVersionString": "1.0",
        "CFBundleVersion": "1", "LSMinimumSystemVersion": "11.0", "NSHighResolutionCapable": True,
        "CFBundleIconFile": "AppIcon", "CFBundleIconName": "AppIcon",
    }
    with open(contents / "Info.plist", "wb") as f:
        plistlib.dump(info, f)
    launcher = LAUNCHER.format(name=name, title=cfg["title"], peel=cfg["peel"], stem=peel.stem, dir=peel.parent,
                               python=Path.home() / "venvs/v/bin/python", script=peel.name)
    tmp = contents / "MacOS" / f".{name}.new"  # replace, never edit in place: a running bash reads its script as it goes
    tmp.write_text(launcher)
    tmp.chmod(0o755)
    tmp.replace(contents / "MacOS" / name)

    for old in [*(contents / "Resources").glob("*.icns"), contents / "Resources" / "Assets.car"]:
        old.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory() as t:
        subprocess.run(["xcrun", "actool", str(pkg), "--compile", str(contents / "Resources"), "--app-icon", "AppIcon",
                        "--platform", "macosx", "--target-device", "mac", "--minimum-deployment-target", "11.0",
                        "--enable-on-demand-resources", "NO", "--development-region", "en", "--include-all-app-icons",
                        "--output-partial-info-plist", f"{t}/partial.plist"], check=True, capture_output=True)
    os.utime(app)  # Finder rereads the icon of a bundle that changed
    subprocess.run([LSREGISTER, "-f", str(app)], check=True)
    return app


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in APPS:
        sys.exit(f"usage: make_app.py {{{'|'.join(APPS)}}}")
    name = sys.argv[1]
    cfg = APPS[name]
    out = OUT / name
    out.mkdir(parents=True, exist_ok=True)
    render(REPO / cfg["peel"], out / "render.png", *cfg["view"])
    art(out / "render.png", out / "art.png", cfg["label"])
    icon_package(out / "art.png", out / "AppIcon.icon", *cfg["fill"])
    print(f"wrote {out / 'render.png'}, {out / 'art.png'}, {out / 'AppIcon.icon'}")
    print(f"built {build_app(name, cfg, out / 'AppIcon.icon')}")


if __name__ == "__main__":
    main()
