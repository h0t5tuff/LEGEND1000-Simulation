#!/usr/bin/env python
"""Build LEGEND.app, one double-click app for the repo's three viewers, and install it in /Applications: l200 and
l1000 (the LEGEND-200 and LEGEND-1000 geometries, peeled layer by layer) and rt (the LEGEND-1000 re-entrance tube
designs). It offers the three when it opens.

    ~/venvs/v/bin/python app/make_app.py

It builds in a temporary folder: each viewer's data (l200 and l1000 the autopeel's meshes, layers and colours; rt
the design scan of output/designs.json and the tube's mesh from the GDML) with its picture for the choices; the icon
(the l1000 picture above the label, an Icon Composer package compiled by actool); legend.swift + autopeel.swift +
rt.swift compiled for Apple silicon and Intel; the bundle, signed ad hoc. It installs it as /Applications/LEGEND.app
and keeps one copy to share, app/output/LEGEND.zip.
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

REPO = Path(__file__).resolve().parents[1]  # app/make_app.py -> the repository
HERE = Path(__file__).resolve().parent
NAME = "LEGEND"
OUT = HERE / "output"
VIEWS = {  # the source of each viewer's geometry, its window title, its picture's camera (theta, phi from +z, deg; zoom[; roll])
    "l200": {"peel": "geom/l200geom/l200-autopeel.py", "title": "LEGEND-200", "view": (62, 25, 1.0)},
    "l1000": {"peel": "geom/l1000geom/l1000-autopeel.py", "title": "LEGEND-1000", "view": (40, 30, 1.15)},
    "rt": {"peel": "geom/l1000geom/l1000-autopeel.py", "title": "LEGEND-1000 re-entrance tube designs", "view": (80, 20, 1.0, -40)},
}
ICON = {"picture": "l1000", "label": "LEGEND", "fill": ("#1d6b6f", "#08262a")}  # background gradient (top, bottom)
SHOW = ("string hardware", "ultem", "PEN")  # autopeel layers drawn with the HPGe detectors in the l200 and l1000 pictures
DESIGNS = REPO / "output" / "designs.json"  # ana/background.C's design scan, which rt shows
STOPS = ((0.9, (38, 70, 150)), (1.0, (33, 120, 190)), (1.1, (40, 170, 150)), (1.2, (140, 200, 80)), (1.3, (250, 220, 60)))  # its colours
SOURCES = ("legend.swift", "autopeel.swift", "rt.swift")
FONT = "/System/Library/Fonts/SFNS.ttf"
MACOS = "12.0"  # oldest macOS the app runs on
LSREGISTER = "/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister"


def load_peel(peel: Path):
    """The autopeel script as a module: its LAYERS, colours and cached geometry, without opening its window."""
    sys.dont_write_bytecode = True  # no __pycache__ left in the geometry folders
    spec = importlib.util.spec_from_file_location("peel", peel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def peel_actors(m) -> list:
    """The HPGe core plus the SHOW layers of an autopeel."""
    layers, _ = m.assemble(m.load_geometry(), m.read_colours())
    return [actor for j, lay in enumerate(layers) if j == len(layers) - 1 or any(k in lay["label"] for k in SHOW)
            for actor, _ in lay["actors"]]


def render(actors: list, png: Path, theta: float, phi: float, zoom: float, roll: float = 0) -> None:
    """The actors, offscreen, on a transparent background."""
    from vtkmodules.vtkIOImage import vtkPNGWriter
    from vtkmodules.vtkRenderingCore import vtkRenderer, vtkRenderWindow, vtkWindowToImageFilter

    ren = vtkRenderer()
    ren.SetBackgroundAlpha(0.0)
    for actor in actors:
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
    cam.Roll(roll)
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


def tube_mesh(m) -> dict:
    """The re-entrant tube, meshed from the GDML's solids finer than the autopeel's (128 slices around, not 16): its SS
    and OFHC Cu layers, and its EFCu solid below the Cu seam (above it the layers are the wall, and the solid's top
    face is no lid). Points, normals [mm, z up] and triangles, all three sections as one mesh."""
    import pyg4ometry
    from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkFiltersCore import vtkClipPolyData, vtkTriangleFilter

    def fine(s):
        if hasattr(s, "nslice"):
            s.nslice = 128
        for b in ("obj1", "obj2"):
            if hasattr(s, b):
                fine(getattr(s, b))

    reg, g = pyg4ometry.gdml.Reader(str(m.GDML)).getRegistry(), m.load_geometry()
    names = [str(n) for n in g["names"]]
    parts = []
    for lv in m.TUBE:  # top to bottom: SS, OFHC Cu, EFCu
        s = reg.logicalVolumeDict[lv].solid
        fine(s)
        v, polys, _ = s.mesh().toVerticesAndPolygons()
        p, n, t = m._surface(np.asarray(v, dtype=float), polys)
        i = names.index(lv)
        assert g["n_pl"][i] == 1, f"{lv} is placed {g['n_pl'][i]} times"
        R, T = g["rot"][g["o_pl"][i]], g["tra"][g["o_pl"][i]]  # its placement in the world
        parts.append([p @ R.T + T, n @ np.linalg.inv(R), t])

    seam = parts[1][0][:, 2].min()  # the Cu seam: the EFCu solid is the wall below it
    pd = vtkPolyData()
    pd.SetPoints(m._points(parts[2][0]))
    cells = vtkCellArray()
    cells.SetData(3, numpy_to_vtkIdTypeArray(parts[2][2].astype(np.int64).ravel(), deep=True))
    pd.SetPolys(cells)
    normals = numpy_to_vtk(np.ascontiguousarray(parts[2][1], dtype=np.float32), deep=True)
    normals.SetName("Normals")
    pd.GetPointData().SetNormals(normals)
    clip = vtkClipPolyData()
    clip.SetInputData(pd)
    clip.SetClipFunction(m.plane((0, 0, seam), (0, 0, -1)))  # keeps z < seam
    tri = vtkTriangleFilter()
    tri.SetInputConnection(clip.GetOutputPort())
    tri.Update()
    out = tri.GetOutput()
    n = vtk_to_numpy(out.GetPointData().GetNormals())
    parts[2] = [vtk_to_numpy(out.GetPoints().GetData()), n / np.linalg.norm(n, axis=1, keepdims=True),
                vtk_to_numpy(out.GetPolys().GetConnectivityArray()).reshape(-1, 3)]

    offset = np.cumsum([0] + [len(p) for p, _, _ in parts])
    return {"pts": np.concatenate([p for p, _, _ in parts]).astype("<f4"), "nrm": np.concatenate([n for _, n, _ in parts]).astype("<f4"),
            "tri": np.concatenate([t + o for (_, _, t), o in zip(parts, offset)]).astype("<u4")}


def tube_actors(tube: dict) -> list:
    """The tube for the icon, painted along its length with the designs' colours: yellow at the top, blue at the bottom."""
    from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper

    p = tube["pts"]
    x = 0.9 + 0.4 * (p[:, 2] - p[:, 2].min()) / np.ptp(p[:, 2])
    xs, cs = np.array([a for a, _ in STOPS]), np.array([c for _, c in STOPS], dtype=float)
    rgb = np.stack([np.interp(x, xs, cs[:, k]) for k in range(3)], axis=1).astype(np.uint8)
    points = vtkPoints()
    points.SetData(numpy_to_vtk(np.ascontiguousarray(p), deep=True))
    pd = vtkPolyData()
    pd.SetPoints(points)
    cells = vtkCellArray()
    cells.SetData(3, numpy_to_vtkIdTypeArray(tube["tri"].astype(np.int64).ravel(), deep=True))
    pd.SetPolys(cells)
    pd.GetPointData().SetNormals(numpy_to_vtk(np.ascontiguousarray(tube["nrm"]), deep=True))
    colours = numpy_to_vtk(rgb, deep=True)
    colours.SetName("rgb")
    pd.GetPointData().SetScalars(colours)
    mapper = vtkPolyDataMapper()
    mapper.SetInputData(pd)
    mapper.SetColorModeToDirectScalars()
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetDiffuse(0.8)
    return [actor]


def outer_radius(rz: tuple, z: np.ndarray) -> np.ndarray:
    """[mm] a closed (r, z) outline's largest radius at each height z (from its crossings, as rz_volume pairs them)."""
    r, zz = rz
    j = np.roll(np.arange(len(r)), 1)
    out = []
    for zm in np.clip(z, zz.min() + 0.01, zz.max() - 0.01):
        cross = (zz > zm) != (zz[j] > zm)
        out.append((r[cross] + (zm - zz[cross]) * (r[j][cross] - r[cross]) / (zz[j][cross] - zz[cross])).max())
    return np.array(out)


def export_designs(m, tube: dict, title: str, out: Path) -> tuple[int, int]:
    """For rt.swift: designs.json (the design scan of output/designs.json, the tube's top and outer radius) and
    tube.bin (the tube's points and normals, float32 x3, then its triangles, uint32 x3)."""
    doc = json.loads(DESIGNS.read_text())
    z = tube["pts"][:, 2]
    if abs(np.ptp(z) - 1000 * doc["depthTot"]) > 1:
        sys.exit(f"the GDML's tube is {np.ptp(z):.0f} mm long, the design scan's {1000 * doc['depthTot']:.0f} mm: rerun ana/background.C")
    step = 10.0  # mm
    depth = np.arange(0, 1000 * doc["depthTot"] + step, step)
    doc |= {"title": title, "zTop": float(z.max()), "rStep": step,
            "radius": [round(float(x), 1) for x in outer_radius(m.gdml_parts()[0][m.TUBE[2]], z.max() - depth)],
            "nPts": len(tube["pts"]), "nTri": len(tube["tri"])}
    with open(out / "tube.bin", "wb") as f:
        for k in ("pts", "nrm", "tri"):
            f.write(tube[k].tobytes())
    (out / "designs.json").write_text(json.dumps(doc, separators=(",", ":")))
    return len(tube["tri"]), (out / "tube.bin").stat().st_size + (out / "designs.json").stat().st_size


def compile_app(exe: Path) -> None:
    """The SOURCES as one universal binary (Apple silicon and Intel)."""
    with tempfile.TemporaryDirectory() as t:
        for arch in ("arm64", "x86_64"):
            subprocess.run(["xcrun", "swiftc", "-O", "-parse-as-library", "-swift-version", "5", "-target",
                            f"{arch}-apple-macos{MACOS}", *[str(HERE / f) for f in SOURCES], "-o", f"{t}/{arch}"], check=True)
        subprocess.run(["lipo", "-create", f"{t}/arm64", f"{t}/x86_64", "-output", str(exe)], check=True)


def build_bundle(app: Path, pkg: Path) -> None:
    """Info.plist, the program and the icon around the data already in Contents/Resources; signed ad hoc."""
    contents = app / "Contents"
    (contents / "MacOS").mkdir(parents=True)
    info = {
        "CFBundleExecutable": NAME, "CFBundleIdentifier": f"local.legend.{NAME}", "CFBundleName": NAME,
        "CFBundleDisplayName": NAME, "CFBundlePackageType": "APPL", "CFBundleShortVersionString": "1.0",
        "CFBundleVersion": "1", "LSMinimumSystemVersion": MACOS, "NSHighResolutionCapable": True,
        "CFBundleIconFile": "AppIcon", "CFBundleIconName": "AppIcon", "NSPrincipalClass": "NSApplication",
        "LSApplicationCategoryType": "public.app-category.education",
    }
    compile_app(contents / "MacOS" / NAME)
    with open(contents / "Info.plist", "wb") as f:
        plistlib.dump(info, f)
    with tempfile.TemporaryDirectory() as t:
        subprocess.run(["xcrun", "actool", str(pkg), "--compile", str(contents / "Resources"), "--app-icon", "AppIcon",
                        "--platform", "macosx", "--target-device", "mac", "--minimum-deployment-target", "11.0",
                        "--enable-on-demand-resources", "NO", "--development-region", "en", "--include-all-app-icons",
                        "--output-partial-info-plist", f"{t}/partial.plist"], check=True, capture_output=True)
    # ad-hoc signature: required to run on Apple silicon, not an identity other Macs trust
    subprocess.run(["xattr", "-cr", str(app)], check=True)  # codesign refuses files with extended attributes
    sign = subprocess.run(["codesign", "--force", "--sign", "-", str(app)], capture_output=True, text=True)
    if sign.returncode:
        sys.exit(f"codesign: {sign.stderr.strip()}")


def main():
    if sys.argv[1:]:
        sys.exit("usage: make_app.py   (no arguments: it builds LEGEND.app with l200, l1000 and rt)")
    if not DESIGNS.exists():
        sys.exit(f"no {DESIGNS}: run ana/background.C first")
    peels = {}
    with tempfile.TemporaryDirectory() as tmp:  # out of iCloud, which puts attributes back that codesign refuses
        t = Path(tmp)
        built, res = t / f"{NAME}.app", t / f"{NAME}.app" / "Contents" / "Resources"
        for key, v in VIEWS.items():
            if v["peel"] not in peels:
                peels[v["peel"]] = load_peel(REPO / v["peel"])
            m = peels[v["peel"]]
            (res / key).mkdir(parents=True)
            if key == "rt":
                tube = tube_mesh(m)
                actors = tube_actors(tube)
                n_tri, n_bytes = export_designs(m, tube, v["title"], res / key)
                print(f"rt: the tube in {n_tri / 1e3:.0f} k triangles; with the design scan {n_bytes / 1e6:.1f} MB")
            else:
                actors = peel_actors(m)
                n_layers, n_bytes = export_geometry(m, v["title"], res / key)
                print(f"{key}: {n_layers} layers, {n_bytes / 1e6:.1f} MB")
            render(actors, res / key / "render.png", *v["view"])  # its picture in the choices
        art(res / ICON["picture"] / "render.png", t / "art.png", ICON["label"])
        icon_package(t / "art.png", t / "AppIcon.icon", *ICON["fill"])
        build_bundle(built, t / "AppIcon.icon")
        OUT.mkdir(exist_ok=True)
        (OUT / f"{NAME}.zip").unlink(missing_ok=True)
        subprocess.run(["ditto", "-c", "-k", "--keepParent", str(built), str(OUT / f"{NAME}.zip")], check=True)
        app = Path("/Applications") / f"{NAME}.app"
        shutil.rmtree(app, ignore_errors=True)  # a running copy keeps its open files; the new one starts next time
        subprocess.run(["ditto", str(built), str(app)], check=True)
    subprocess.run(["xattr", "-cr", str(app)], check=True)  # the copy picks up Finder info, which fails a strict check
    os.utime(app)  # Finder rereads the icon of a bundle that changed
    subprocess.run([LSREGISTER, "-f", str(app)], check=True)
    print(f"installed {app}; to share: {OUT / f'{NAME}.zip'}")


if __name__ == "__main__":
    main()
