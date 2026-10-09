#!/usr/bin/env python
"""Build a double-click app for one of the repo's autopeels, with an icon of its geometry, into /Applications.

    ~/venvs/v/bin/python geom/apps/make_app.py l200              (or l1000; see README.md)
    ~/venvs/v/bin/python geom/apps/make_app.py l200 --launcher   the old kind: runs the repo's autopeel

1. render the detector array offscreen from the autopeel's geometry   -> output/<name>/render.png
2. put the label under it, on a transparent 1024 px canvas            -> output/<name>/art.png
3. wrap it in an Icon Composer package with the background gradient   -> output/<name>/AppIcon.icon
4. standalone (default): export the autopeel's meshes, layers and colours (geometry.json, geometry.bin), compile
   autopeel.swift for Apple silicon and Intel, bundle, sign ad hoc -> output/<name>/<name>.app and <name>.zip;
   --launcher: a bundle whose bash script runs the autopeel with ~/venvs/v
5. install it as /Applications/<name>.app (replacing the one there), its icon compiled by Xcode's actool
"""

import importlib.util
import json
import os
import plistlib
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

REPO = Path(__file__).resolve().parents[2]  # geom/apps/make_app.py -> the repository
HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
APPS = {  # autopeel, title, label, camera (theta, phi from +z, deg; zoom), background gradient (top, bottom)
    "l200": {"peel": "LEGEND200/l200-autopeel.py", "title": "LEGEND-200", "label": "L200", "view": (62, 25, 1.0), "fill": ("#2b4f86", "#0d1830")},
    "l1000": {"peel": "geom/l1000-autopeel.py", "title": "LEGEND-1000", "label": "L1000", "view": (40, 30, 1.15), "fill": ("#1d6b6f", "#08262a")},
}
SHOW = ("string hardware", "ultem", "PEN")  # autopeel layers drawn with the HPGe detectors in the icon
FONT = "/System/Library/Fonts/SFNS.ttf"
MACOS = "12.0"  # oldest macOS the standalone app runs on
LSREGISTER = "/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister"


def load_peel(peel: Path):
    """The autopeel script as a module: its LAYERS, colours and cached geometry, without opening its window."""
    spec = importlib.util.spec_from_file_location("peel", peel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def render(m, png: Path, theta: float, phi: float, zoom: float) -> None:
    """The HPGe core plus the SHOW layers of an autopeel, offscreen, on a transparent background."""
    from vtkmodules.vtkIOImage import vtkPNGWriter
    from vtkmodules.vtkRenderingCore import vtkRenderer, vtkRenderWindow, vtkWindowToImageFilter

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


def export_geometry(m, title: str, out: Path) -> tuple[int, int]:
    """The autopeel's meshes (one per logical volume), their placements, its layers and colours, for autopeel.swift:
    geometry.json (layers -> groups of one colour and hierarchy depth -> parts) and geometry.bin (the arrays)."""
    g, colours = m.load_geometry(), m.read_colours()
    labels = [lay["label"] for lay in m.assemble(g, colours)[0]]  # the autopeel's own layer names, its core's last
    names = [str(n) for n in g["names"]]
    regexes = [re.compile(rx) for _, rx in m.LAYERS]
    groups = [{} for _ in labels]  # the same sorting as the autopeel's assemble()
    for i, name in enumerate(names):
        rgba = tuple(float(c) for c in colours.get(name, (0.7, 0.7, 0.7, 1.0)))
        if rgba[3] == 0:  # the world box
            continue
        j = next((j for j, rx in enumerate(regexes) if rx.fullmatch(name)), len(m.LAYERS))
        groups[j].setdefault((rgba, int(g["depth"][i])), []).append(i)

    o = {k: np.concatenate([[0], np.cumsum(g[f"n_{k}"])]) for k in ("pts", "tri", "pl")}
    arrays = {"pts": g["pts"].astype("<f4"), "nrm": g["nrm"].astype("<f4"), "tri": g["tri"].astype("<u4"),
              "rot": g["rot"].reshape(-1, 9).astype("<f4"), "tra": g["tra"].astype("<f4")}
    index, offset = {}, 0
    with open(out / "geometry.bin", "wb") as f:
        for k, a in arrays.items():
            index[k] = [offset, int(a.size)]
            f.write(a.tobytes())
            offset += a.nbytes
    doc = {
        "title": title, "start_theta_phi": list(m.START_THETA_PHI), "dwell": m.DWELL, "fade": m.FADE,
        "layers": [{"label": lab, "groups": [{"rgba": list(rgba), "depth": d, "parts": parts}
                                              for (rgba, d), parts in grp.items()]} for lab, grp in zip(labels, groups)],
        "parts": [{"name": n, "pts": [int(o["pts"][i]), int(g["n_pts"][i])], "tri": [int(o["tri"][i]), int(g["n_tri"][i])],
                   "pl": [int(o["pl"][i]), int(g["n_pl"][i])]} for i, n in enumerate(names)],
        "arrays": index,
    }
    (out / "geometry.json").write_text(json.dumps(doc, separators=(",", ":")))
    return len(labels), offset


def compile_viewer(exe: Path) -> None:
    """autopeel.swift as one universal binary (Apple silicon and Intel)."""
    with tempfile.TemporaryDirectory() as t:
        for arch in ("arm64", "x86_64"):
            subprocess.run(["xcrun", "swiftc", "-O", "-parse-as-library", "-swift-version", "5", "-target",
                            f"{arch}-apple-macos{MACOS}", str(HERE / "autopeel.swift"), "-o", f"{t}/{arch}"], check=True)
        subprocess.run(["lipo", "-create", f"{t}/arm64", f"{t}/x86_64", "-output", str(exe)], check=True)


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


def build_bundle(name: str, cfg: dict, m, pkg: Path, app: Path, standalone: bool) -> None:
    """A fresh <name>.app at `app`: Info.plist, the program (viewer + geometry, or the launcher), the icon."""
    shutil.rmtree(app, ignore_errors=True)
    contents, peel = app / "Contents", REPO / cfg["peel"]
    (contents / "MacOS").mkdir(parents=True)
    (contents / "Resources").mkdir()
    info = {
        "CFBundleExecutable": name, "CFBundleIdentifier": f"local.legend.{name}", "CFBundleName": name,
        "CFBundleDisplayName": name, "CFBundlePackageType": "APPL", "CFBundleShortVersionString": "1.0",
        "CFBundleVersion": "1", "LSMinimumSystemVersion": MACOS if standalone else "11.0", "NSHighResolutionCapable": True,
        "CFBundleIconFile": "AppIcon", "CFBundleIconName": "AppIcon",
    }
    if standalone:
        info |= {"NSPrincipalClass": "NSApplication", "LSApplicationCategoryType": "public.app-category.education"}
        compile_viewer(contents / "MacOS" / name)
        n_layers, n_bytes = export_geometry(m, cfg["title"], contents / "Resources")
        print(f"geometry: {n_layers} layers, {n_bytes / 1e6:.1f} MB")
    else:
        (contents / "MacOS" / name).write_text(LAUNCHER.format(
            name=name, title=cfg["title"], peel=cfg["peel"], stem=peel.stem, dir=peel.parent,
            python=Path.home() / "venvs/v/bin/python", script=peel.name))
        (contents / "MacOS" / name).chmod(0o755)
    with open(contents / "Info.plist", "wb") as f:
        plistlib.dump(info, f)
    with tempfile.TemporaryDirectory() as t:
        subprocess.run(["xcrun", "actool", str(pkg), "--compile", str(contents / "Resources"), "--app-icon", "AppIcon",
                        "--platform", "macosx", "--target-device", "mac", "--minimum-deployment-target", "11.0",
                        "--enable-on-demand-resources", "NO", "--development-region", "en", "--include-all-app-icons",
                        "--output-partial-info-plist", f"{t}/partial.plist"], check=True, capture_output=True)
    if standalone:  # ad-hoc signature: required to run on Apple silicon, not an identity other Macs trust
        subprocess.run(["xattr", "-cr", str(app)], check=True)  # codesign refuses files with extended attributes
        subprocess.run(["codesign", "--force", "--sign", "-", str(app)], check=True, capture_output=True)


def main():
    args = sys.argv[1:]
    standalone = "--launcher" not in args
    names = [a for a in args if a != "--launcher"]
    if len(names) != 1 or names[0] not in APPS:
        sys.exit(f"usage: make_app.py {{{'|'.join(APPS)}}} [--launcher]")
    name = names[0]
    cfg = APPS[name]
    out = OUT / name
    out.mkdir(parents=True, exist_ok=True)
    m = load_peel(REPO / cfg["peel"])
    render(m, out / "render.png", *cfg["view"])
    art(out / "render.png", out / "art.png", cfg["label"])
    icon_package(out / "art.png", out / "AppIcon.icon", *cfg["fill"])
    print(f"wrote {out / 'render.png'}, {out / 'art.png'}, {out / 'AppIcon.icon'}")

    built = out / f"{name}.app"
    build_bundle(name, cfg, m, out / "AppIcon.icon", built, standalone)
    if standalone:
        (out / f"{name}.zip").unlink(missing_ok=True)
        subprocess.run(["ditto", "-c", "-k", "--keepParent", str(built), str(out / f"{name}.zip")], check=True)
        print(f"wrote {built} and {out / f'{name}.zip'}")
    app = Path("/Applications") / f"{name}.app"
    shutil.rmtree(app, ignore_errors=True)  # a running copy keeps its open files; the new one starts next time
    subprocess.run(["ditto", str(built), str(app)], check=True)
    if standalone:
        subprocess.run(["xattr", "-cr", str(app)], check=True)  # the copy picks up Finder info, which fails a strict check
    else:
        shutil.rmtree(built)
    os.utime(app)  # Finder rereads the icon of a bundle that changed
    subprocess.run([LSREGISTER, "-f", str(app)], check=True)
    print(f"installed {app} ({'standalone' if standalone else 'launcher'})")


if __name__ == "__main__":
    main()
