#!/usr/bin/env python
"""Peel the LEGEND-1000 geometry layer by layer on a timer while you rotate it (see README.md).

    ~/venvs/v/bin/python l1000-autopeel.py
    ~/venvs/v/bin/python l1000-autopeel.py --holders ../output/holders_geometry.png    (the detector holders alone, no window)
    ~/venvs/v/bin/python l1000-autopeel.py --tube ../output/tube_geometry.png          (the re-entrant tube, the decay source)
    ~/venvs/v/bin/python l1000-autopeel.py --uglar ../output/uglar_geometry.png        (the underground argon, the veto)
    ~/venvs/v/bin/python l1000-autopeel.py --hpge ../output/hpge_geometry.png          (the HPGe detectors)
"""

import itertools
import re
import sys
import time
from pathlib import Path

import numpy as np
import vtkmodules.vtkRenderingFreeType  # noqa: F401  (text)
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401  (the OpenGL render window)
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkMapper,
    vtkPolyDataMapper,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkTextActor,
)

HERE = Path(__file__).resolve().parent
GDML = HERE / "l1000.gdml"
COLOURS = HERE / "l1000-colors.mac"
CACHE = HERE / "autopeel-cache.npz"
DWELL = 3.0  # s per layer at start (up/down change it)
FADE = 0.5  # s a layer takes to fade out or back in
START_THETA_PHI = (70, 20)  # deg, z up: the start view of l1000-vis.mac

# Peel order, outside in. Each LV goes to the first layer whose regex matches its whole name;
# whatever no layer matches (the HPGe detectors V00101Z ...) is left at the end.
LAYERS = [
    ("outer cryostat: outer steel tank, skirt and foot", r"cryostat_(outer|skirt|foot)_steel_316L"),
    ("vacuum insulation between the two steel tanks", r"cryostat_insulation_vacuum"),
    ("inner cryostat", r"cryostat_inner_steel_316L"),
    ("atmospheric LAr (outer bath)", r"liquid_argon_atmospheric"),
    ("neutron moderator", r"neutron_moderator_pmma"),
    ("outer reflector: TPB coating", r"atmospheric_wlsr_tpb"),
    ("outer reflector: tetratex lining", r"atmospheric_wlsr_tetratex"),
    ("re-entrant tube: EFCu (z < 920 mm), OFHC Cu (920-2925), SS 316L (> 2925)", r"reentrance_tube_.*"),
    ("underground LAr (inner active bath)", r"liquid_argon_underground"),
    ("inner reflector: TPB coating", r"underground_wlsr_tpb"),
    ("inner reflector: tetratex lining", r"underground_wlsr_tetratex"),
    ("WLS fibres: outer and inner cladding, TPB coating, PS core", r"fiber_.*"),
    ("SiPM mounts (copper rings, top and bottom)", r"larinstr_.*_sipm_wrap_copper_.*"),
    ("SiPMs", r"larinstr_.*_sipm_silicon_.*"),
    ("string hardware: hangers, rods, tristars, weldments, cable caps, cables", r"hpge_string_.*|hpge_cable_.*"),
    ("ultem insulators and signal/HV clamps", r"hpge_assembly_(clamp_.*|insulator)_ultem.*"),
    ("PEN holder plates", r"hpge_assembly_plate_pen_.*"),
    ("signal ASICs (1 mm silica chips under each detector: look from below)", r"hpge_assembly_asic"),
]


# ── geometry: GDML -> one surface per LV + its placements (cached) ──────────────────────

def _points(a):
    pts = vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(a, dtype=np.float32), deep=True))
    return pts


def _surface(verts, polys):
    """One LV's CSG polygons -> triangles with smooth normals, split at edges sharper than 30 deg."""
    from vtkmodules.vtkFiltersCore import vtkCleanPolyData, vtkPolyDataNormals, vtkTriangleFilter

    offsets = np.zeros(len(polys) + 1, dtype=np.int64)
    np.cumsum([len(p) for p in polys], out=offsets[1:])
    conn = np.fromiter(itertools.chain.from_iterable(polys), dtype=np.int64, count=offsets[-1])
    cells = vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(offsets, deep=True), numpy_to_vtkIdTypeArray(conn, deep=True))
    pd = vtkPolyData()
    pd.SetPoints(_points(verts))
    pd.SetPolys(cells)

    clean = vtkCleanPolyData()  # CSG output repeats each vertex per polygon: merge them
    clean.SetInputData(pd)
    tri = vtkTriangleFilter()
    tri.SetInputConnection(clean.GetOutputPort())
    tri.PassLinesOff()
    tri.PassVertsOff()
    nrm = vtkPolyDataNormals()
    nrm.SetInputConnection(tri.GetOutputPort())
    nrm.SetFeatureAngle(30)
    nrm.SplittingOn()
    nrm.ConsistencyOff()  # keep the CSG winding (outward)
    nrm.ComputeCellNormalsOff()
    nrm.Update()
    out = nrm.GetOutput()
    if out.GetNumberOfPolys() == 0:
        return None
    return (vtk_to_numpy(out.GetPoints().GetData()).astype(np.float32),
            vtk_to_numpy(out.GetPointData().GetNormals()).astype(np.float32),
            vtk_to_numpy(out.GetPolys().GetConnectivityArray()).reshape(-1, 3).astype(np.int32))


def build_cache(key):
    print(f"meshing {GDML.name} (first run, or the GDML changed): ~25 s", flush=True)
    from pyg4ometry import gdml
    from pyg4ometry.visualisation.ViewerBase import ViewerBase

    reg = gdml.Reader(str(GDML)).getRegistry()  # meshes every LV on the way
    walk = ViewerBase()
    walk.addLogicalVolume(reg.worldVolume)  # world transform of every placement, per LV

    depth = {reg.worldVolume.name: 0}  # hierarchy level, for the coincident-face offset

    def tree(lv):
        for pv in lv.daughterVolumes:
            if pv.logicalVolume.name not in depth:
                depth[pv.logicalVolume.name] = depth[lv.name] + 1
                tree(pv.logicalVolume)

    tree(reg.worldVolume)

    out = {k: [] for k in ("names", "depth", "pts", "nrm", "tri", "rot", "tra", "n_pts", "n_tri", "n_pl")}
    for name, csg in walk.localmeshes.items():
        if name not in reg.logicalVolumeDict:  # overlap / extrusion helper meshes
            continue
        verts, polys, _ = csg.toVerticesAndPolygons()
        surf = _surface(np.asarray(verts, dtype=float), polys)
        if surf is None:
            continue
        places = walk.instancePlacements[name]
        out["names"].append(name)
        out["depth"].append(depth.get(name, 0))
        for k, a in zip(("pts", "nrm", "tri"), surf):
            out[k].append(a)
        out["rot"].append(np.array([np.asarray(p["transformation"], dtype=float) for p in places]))
        out["tra"].append(np.array([np.asarray(p["translation"], dtype=float).ravel() for p in places]))
        out["n_pts"].append(len(surf[0]))
        out["n_tri"].append(len(surf[2]))
        out["n_pl"].append(len(places))

    g = {k: (np.concatenate(v) if k in ("pts", "nrm", "tri", "rot", "tra") else np.array(v))
         for k, v in out.items()}
    g["key"] = np.array(key)
    CACHE.parent.mkdir(exist_ok=True)
    tmp = CACHE.with_suffix(".tmp")
    with open(tmp, "wb") as f:  # write then rename, so an interrupted build leaves no bad cache
        np.savez(f, **g)
    tmp.replace(CACHE)
    return g


def load_geometry():
    st = GDML.stat()
    key = f"{GDML.resolve()}|{st.st_size}|{st.st_mtime_ns}"
    g = None
    if CACHE.exists():
        with np.load(CACHE) as z:
            if str(z["key"]) == key:
                g = {k: z[k] for k in z.files}
    if g is None:
        g = build_cache(key)
    for k in ("pts", "tri", "pl"):
        g[f"o_{k}"] = np.concatenate([[0], np.cumsum(g[f"n_{k}"])])
    return g


def placed(g, i, keep=None):
    """LV i's surface at all of its placements (or those `keep` selects), in world coordinates."""
    p = g["pts"][g["o_pts"][i]:g["o_pts"][i + 1]]
    n = g["nrm"][g["o_pts"][i]:g["o_pts"][i + 1]]
    t = g["tri"][g["o_tri"][i]:g["o_tri"][i + 1]]
    rot = g["rot"][g["o_pl"][i]:g["o_pl"][i + 1]]
    tra = g["tra"][g["o_pl"][i]:g["o_pl"][i + 1]]
    if keep is not None:
        rot, tra = rot[keep], tra[keep]
    pts = np.einsum("kij,nj->kni", rot, p) + tra[:, None, :]
    nrm = np.einsum("kij,nj->kni", np.linalg.inv(rot).transpose(0, 2, 1), n)  # normals: inverse transpose
    nrm /= np.maximum(np.linalg.norm(nrm, axis=2, keepdims=True), 1e-12)
    tri = t[None] + (np.arange(len(rot)) * len(p))[:, None, None]
    return pts.reshape(-1, 3), nrm.reshape(-1, 3), tri.reshape(-1, 3)


# ── scene: one actor per (layer, colour, depth) ─────────────────────────────────────────

def read_colours():
    rgba = {}
    for line in COLOURS.read_text().splitlines():
        f = line.split()
        if len(f) == 7 and f[0] == "/vis/geometry/set/colour":
            rgba[f[1]] = tuple(float(x) for x in f[3:])
    return rgba


def make_actor(parts, rgba, depth, clip=None, cap=False):
    """One actor of the given parts. clip: a vtkPlane that cuts away the half its normal points from;
    cap: close the cut faces (for a closed solid, e.g. a detector), else leave the cut open (thin walls)."""
    n_pts = np.cumsum([0] + [len(p) for p, _, _ in parts])
    pts = np.concatenate([p for p, _, _ in parts])
    nrm = np.concatenate([n for _, n, _ in parts])
    tri = np.concatenate([t + o for (_, _, t), o in zip(parts, n_pts)])
    pd = vtkPolyData()
    pd.SetPoints(_points(pts))
    cells = vtkCellArray()
    cells.SetData(3, numpy_to_vtkIdTypeArray(tri.astype(np.int64).ravel(), deep=True))
    pd.SetPolys(cells)
    normals = numpy_to_vtk(np.ascontiguousarray(nrm, dtype=np.float32), deep=True)
    normals.SetName("Normals")
    pd.GetPointData().SetNormals(normals)

    mapper = vtkPolyDataMapper()
    if clip is not None and cap:
        from vtkmodules.vtkCommonDataModel import vtkPlaneCollection
        from vtkmodules.vtkFiltersCore import vtkPolyDataNormals
        from vtkmodules.vtkFiltersGeneral import vtkClipClosedSurface

        planes = vtkPlaneCollection()
        planes.AddItem(clip)
        cut = vtkClipClosedSurface()
        cut.SetInputData(pd)
        cut.SetClippingPlanes(planes)
        nrm_f = vtkPolyDataNormals()  # the caps need normals of their own
        nrm_f.SetInputConnection(cut.GetOutputPort())
        nrm_f.SetFeatureAngle(30)
        nrm_f.SplittingOn()
        nrm_f.Update()
        mapper.SetInputConnection(nrm_f.GetOutputPort())
    else:
        mapper.SetInputData(pd)
        if clip is not None:
            mapper.AddClippingPlane(clip)
    mapper.ScalarVisibilityOff()
    mapper.StaticOn()  # geometry never changes: upload once
    # a daughter that shares a face with its mother (e.g. the OFHC and SS layers on the re-entrant
    # tube) wins the depth test, instead of the two z-fighting
    mapper.SetRelativeCoincidentTopologyPolygonOffsetParameters(0, -2.0 * depth)
    actor = vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetColor(*rgba[:3])
    prop.SetOpacity(rgba[3])
    prop.SetAmbient(0.15)
    prop.SetDiffuse(0.85)
    prop.SetSpecular(0.15)
    prop.SetSpecularPower(25)
    return actor, len(tri)


def assemble(g, colours):
    """Sort the LVs into LAYERS plus the core that is never peeled, and build their actors.
    Returns (layers, triangles); the core is the last layer."""
    names = [str(n) for n in g["names"]]
    missing = [n for n in names if n not in colours]
    if missing:
        print(f"warning: {len(missing)} LVs have no colour in {COLOURS.name} (drawn grey), e.g. {missing[:3]}")
    regexes = [re.compile(rx) for _, rx in LAYERS]
    groups = [{} for _ in range(len(LAYERS) + 1)]  # per layer: (rgba, depth) -> [placed parts]
    counts = [[0, 0] for _ in groups]  # per layer: LVs, placements
    left, n_tri = [], 0
    for i, name in enumerate(names):
        rgba = colours.get(name, (0.7, 0.7, 0.7, 1.0))
        if rgba[3] == 0:  # the world box
            continue
        layer = next((j for j, rx in enumerate(regexes) if rx.fullmatch(name)), len(LAYERS))
        if layer == len(LAYERS):
            left.append(name)
        groups[layer].setdefault((rgba, int(g["depth"][i])), []).append(placed(g, i))
        counts[layer][0] += 1
        counts[layer][1] += int(g["n_pl"][i])

    core = (f"the {len(left)} HPGe detectors" if all(re.fullmatch(r"V\d{5}Z", n) for n in left)
            else f"{len(left)} volumes no layer matches")
    layers = []
    for j, grp in enumerate(groups):
        actors = []
        for (rgba, depth), parts in grp.items():
            actor, nt = make_actor(parts, rgba, depth)
            actors.append((actor, rgba[3]))  # with its own opacity, which fades scale
            n_tri += nt
        label = LAYERS[j][0] if j < len(LAYERS) else core
        layers.append({"label": label, "actors": actors, "lvs": counts[j][0], "placements": counts[j][1]})
    return layers, n_tri


# ── the peel: a timer that only changes visibility and opacity ──────────────────────────

class Peeler:
    def __init__(self, layers, win, ren, text):
        self.layers, self.win, self.ren, self.text = layers, win, ren, text
        self.n = len(layers) - 1  # the last entry is the core, which stays
        self.alpha = [1.0] * self.n  # current and target fade of each layer
        self.target = [1.0] * self.n
        self.peeled, self.dwell, self.playing = 0, DWELL, True
        self.last = time.monotonic()
        self.next_at = self.last + DWELL
        self.dirty = True  # redraw on the next tick

    def tick(self):
        now = time.monotonic()
        dt, self.last = now - self.last, now
        if self.playing and now >= self.next_at:
            self.forward()
            self.dirty = True
        if self.animate(dt) or self.dirty:
            self.dirty = False
            self.show_text()
            self.ren.ResetCameraClippingRange()  # what is visible changed, the camera did not
            self.win.Render()

    def animate(self, dt):
        moved = False
        for i in range(self.n):
            a, t = self.alpha[i], self.target[i]
            if a == t:
                continue
            a = max(t, a - dt / FADE) if t < a else min(t, a + dt / FADE)
            self.alpha[i] = a
            for actor, opacity in self.layers[i]["actors"]:
                actor.SetVisibility(a > 0)
                actor.GetProperty().SetOpacity(a * opacity)
            moved = True
        return moved

    def forward(self):
        if self.peeled < self.n:
            self.target[self.peeled] = 0.0
            self.peeled += 1
        self.playing = self.playing and self.peeled < self.n
        self.next_at = time.monotonic() + self.dwell

    def back(self):
        if self.peeled > 0:
            self.peeled -= 1
            self.target[self.peeled] = 1.0
        self.next_at = time.monotonic() + self.dwell

    def restart(self):
        self.peeled, self.target, self.playing = 0, [1.0] * self.n, True
        self.next_at = time.monotonic() + self.dwell

    def toggle(self):
        if self.peeled == self.n:
            return self.restart()
        self.playing = not self.playing
        self.next_at = time.monotonic() + self.dwell

    def speed(self, factor):
        self.dwell = min(max(self.dwell * factor, 0.3), 60.0)
        self.next_at = min(self.next_at, time.monotonic() + self.dwell)

    def show_text(self):
        state = f"peeling every {self.dwell:.1f} s" if self.playing else "paused"
        if self.peeled == self.n:
            state = "done (space starts again)"
        lines = [f"{self.peeled}/{self.n} layers removed    {state}"]
        if self.peeled:
            lines.append(f"removed:  {self.layers[self.peeled - 1]['label']}")
        nxt = self.layers[self.peeled]
        lines.append(f"{'next' if self.peeled < self.n else 'left'}:  {nxt['label']}")
        self.text.SetInput("\n".join(lines))


# ── window ──────────────────────────────────────────────────────────────────────────────

def start_view(ren, bounds):
    th, ph = np.radians(START_THETA_PHI)
    cam = ren.GetActiveCamera()
    cam.SetFocalPoint(0, 0, 0)
    cam.SetPosition(np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th))
    cam.SetViewUp(0, 0, 1)
    ren.ResetCamera(bounds)  # keeps the direction, frames the whole detector


def text_actor(size, colour, x, y):
    """Text anchored at (x, y) in window fractions; it hangs down from the top half, sits up from the bottom."""
    t = vtkTextActor()
    p = t.GetTextProperty()
    p.SetFontSize(size)
    p.SetColor(*colour)
    p.SetShadow(True)
    p.SetLineSpacing(1.3)
    t.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    t.SetPosition(x, y)
    if y > 0.5:
        p.SetVerticalJustificationToTop()
    else:
        p.SetVerticalJustificationToBottom()
    return t


# ── the detector holders alone, offscreen to a PNG (--holders) ──────────────────────────

HOLDERS = "hpge_string_support_weldment_copper"  # legend-pygeom-l1000: "the copper weldment holding a detector unit to the support rods"
UNIT = "V00101Z"  # the close-up's detector


def holders_png(out):
    """Left: every detector holder alone. Right: one detector unit, its holders solid, the detector and PEN plate ghosted."""
    g = load_geometry()
    names = [str(n) for n in g["names"]]
    h, d = names.index(HOLDERS), names.index(UNIT)
    tra_h = g["tra"][g["o_pl"][h]:g["o_pl"][h + 1]]
    c = g["tra"][g["o_pl"][d]]  # the detector's origin, at its bottom
    near = (np.abs(tra_h[:, 2] - c[2]) < 30) & (np.hypot(*(tra_h[:, :2] - c[:2]).T) < 150)  # mm: this unit's holders
    pen = [i for i, n in enumerate(names) if n.startswith("hpge_assembly_plate_pen")]
    pen_keep = {i: np.linalg.norm(g["tra"][g["o_pl"][i]:g["o_pl"][i + 1]] - c, axis=1) < 30 for i in pen}

    try:
        import yaml
        part = next(p for p in yaml.safe_load((HERE / "l1000-parts.yaml").read_text())["parts"] if p["name"] == HOLDERS)
        what = f"{part['material']}, {part['unit_volume'] * 1e3:.0f} mm3 each, {part['total_mass'] / 1e3:.2f} kg in all"
    except (OSError, StopIteration):
        what = "(no l1000-parts.yaml)"
    print(f"the detector holders: {HOLDERS}, {int(g['n_pl'][h])} placements, 3 per detector ({int(g['n_pl'][h]) // 3} detectors); {what}")
    print(f"close-up: {UNIT} at ({c[0]:.0f}, {c[1]:.0f}, {c[2]:.0f}) mm, its {int(near.sum())} holders "
          f"{c[2] - tra_h[near, 2].max():.1f} mm below the detector's base")

    copper = (0.85, 0.50, 0.25, 1.0)
    write_panels(out, [
        (f"all {int(g['n_pl'][h])} EFCu detector holders, alone", [make_actor([placed(g, h)], copper, 0)[0]], (*START_THETA_PHI, 1.35)),
        (f"one detector unit ({UNIT}): its 3 holders; the detector and its PEN plate ghosted",
         [make_actor([placed(g, h, near)], copper, 0)[0], make_actor([placed(g, d)], (0.75, 0.80, 0.85, 0.25), 0)[0],
          make_actor([placed(g, i, m) for i, m in pen_keep.items() if m.any()], (0.55, 0.75, 0.95, 0.25), 0)[0]],
         (*START_THETA_PHI, 1.0)),
    ])


def write_panels(out, panels):
    """Two panels side by side, offscreen, to a 2400 x 1200 PNG. Each panel: (label, actors, (theta, phi, zoom)
    or (theta, phi, zoom, bounds)); the camera looks from theta/phi (deg, z up) and frames the bounds (default: its actors)."""
    vtkMapper.SetResolveCoincidentTopologyToPolygonOffset()
    win = vtkRenderWindow()
    win.SetOffScreenRendering(1)
    win.SetSize(2400, 1200)
    for k, (label, actors, view) in enumerate(panels):
        ren = vtkRenderer()
        ren.SetViewport(0.5 * k, 0, 0.5 * k + 0.5, 1)
        ren.GradientBackgroundOn()
        ren.SetBackground(0.05, 0.05, 0.07)
        ren.SetBackground2(0.20, 0.22, 0.26)
        ren.SetUseDepthPeeling(True)
        for a in actors:
            ren.AddActor(a)
        th, ph = np.radians(view[:2])
        cam = ren.GetActiveCamera()
        cam.SetFocalPoint(0, 0, 0)
        cam.SetPosition(np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th))
        cam.SetViewUp(0, 0, 1)
        ren.ResetCamera(view[3] if len(view) > 3 else ren.ComputeVisiblePropBounds())
        cam.Zoom(view[2])
        ren.ResetCameraClippingRange()
        t = text_actor(28, (0.95, 0.95, 0.95), 0.02, 0.97)
        t.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()  # this panel's corner, not the window's
        t.SetInput(label)
        ren.AddViewProp(t)
        win.AddRenderer(ren)
    win.Render()

    from vtkmodules.vtkIOImage import vtkPNGWriter
    from vtkmodules.vtkRenderingCore import vtkWindowToImageFilter
    grab = vtkWindowToImageFilter()
    grab.SetInput(win)
    grab.Update()
    png = vtkPNGWriter()
    png.SetFileName(str(out))
    png.SetInputConnection(grab.GetOutputPort())
    png.Write()
    print(f"wrote {out}")


# ── the tube, the UGLAr and the HPGe detectors of sim/run.mac, offscreen to PNGs (--tube, --uglar, --hpge) ──

TUBE = ("reentrance_tube_layer_steel_316L", "reentrance_tube_layer_copper_ofhc", "reentrance_tube_copper")  # top to bottom
UGLAR = "liquid_argon_underground"
GHOST = (0.75, 0.80, 0.85, 0.25)


def parts_info(names):
    """name -> its line of l1000-parts.yaml (material, placements, volume, mass), if the manifest is there."""
    try:
        import yaml
        return {p["name"]: p for p in yaml.safe_load((HERE / "l1000-parts.yaml").read_text())["parts"] if p["name"] in names}
    except OSError:
        return {}


def plane(origin, normal):
    from vtkmodules.vtkCommonDataModel import vtkPlane
    p = vtkPlane()
    p.SetOrigin(*origin)
    p.SetNormal(*normal)
    return p


def facing(phi, origin=(0, 0, 0)):
    """The vertical plane through `origin` facing the camera at azimuth phi: it cuts away the half nearer the camera."""
    return plane(origin, (-np.cos(np.radians(phi)), -np.sin(np.radians(phi)), 0))


def _bounds(pts, zmax=None):
    if zmax is not None:
        pts = pts[pts[:, 2] < zmax]
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    return [lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]]


def half_solid(name, origin, phi, rgba):
    """A GDML polycone (a detector) as a solid half: its exact (r, z) profile swept 180 deg about its axis, both cut
    faces closed, the half towards the camera at azimuth phi left out; placed at origin (the detectors are not rotated)."""
    from vtkmodules.vtkFiltersCore import vtkPolyDataNormals, vtkTriangleFilter
    from vtkmodules.vtkFiltersModeling import vtkRotationalExtrusionFilter

    block = re.search(rf'<genericPolycone name="{name}"[^>]*>(.*?)</genericPolycone>', GDML.read_text(), re.S).group(1)
    rz = [(float(r), float(z)) for r, z in re.findall(r'<rzpoint r="([^"]+)" z="([^"]+)"', block)]
    a0 = np.radians(phi) + np.pi / 2  # the half kept spans phi + 90 .. phi + 270 deg, away from the camera
    pts, cell = vtkPoints(), vtkCellArray()
    cell.InsertNextCell(len(rz))
    for k, (r, z) in enumerate(rz):
        pts.InsertNextPoint(r * np.cos(a0), r * np.sin(a0), z)
        cell.InsertCellPoint(k)
    pd = vtkPolyData()
    pd.SetPoints(pts)
    pd.SetPolys(cell)
    sweep = vtkRotationalExtrusionFilter()
    sweep.SetInputData(pd)
    sweep.SetAngle(180)
    sweep.SetResolution(90)
    sweep.CappingOn()  # the two cut faces: copies of the profile
    tri = vtkTriangleFilter()
    tri.SetInputConnection(sweep.GetOutputPort())
    nrm = vtkPolyDataNormals()
    nrm.SetInputConnection(tri.GetOutputPort())
    nrm.SetFeatureAngle(30)
    nrm.SplittingOn()
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(nrm.GetOutputPort())
    mapper.ScalarVisibilityOff()
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.SetPosition(*origin)
    prop = actor.GetProperty()
    prop.SetColor(*rgba[:3])
    prop.SetAmbient(0.15)
    prop.SetDiffuse(0.85)
    prop.SetSpecular(0.15)
    prop.SetSpecularPower(25)
    return actor


def tube_png(out):
    """Left: the tube's three sections, the decay sources. Right: its bottom cut along the axis, the detectors ghosted inside."""
    g, colours = load_geometry(), read_colours()
    names = [str(n) for n in g["names"]]
    idx = [names.index(n) for n in TUBE]
    dets = [i for i, n in enumerate(names) if re.fullmatch(r"V\d{5}Z", n)]
    info = parts_info(TUBE)
    pts = np.concatenate([placed(g, i)[0] for i in idx])
    seams = [placed(g, i)[0][:, 2].min() for i in idx[:2]]
    for n in TUBE:
        p = info.get(n, {})
        print(f"{n}: {p.get('material', '?')}, {p.get('total_mass', 0) / 1e3:.1f} kg")
    print(f"z {pts[:, 2].min():.0f} .. {pts[:, 2].max():.0f} mm; seams: steel above z {seams[0]:.0f}, OFHC above z {seams[1]:.0f}; "
          f"the detectors at z {min(placed(g, i)[0][:, 2].min() for i in dets):.0f} .. {max(placed(g, i)[0][:, 2].max() for i in dets):.0f}")
    tube = lambda clip=None: [make_actor([placed(g, i)], colours[n], int(g["depth"][i]), clip)[0] for i, n in zip(idx, TUBE)]
    th, ph = START_THETA_PHI
    write_panels(out, [
        ("the re-entrant tube, sim/run.mac's decay source: steel, OFHC Cu, EFCu", tube(), (th, ph, 0.9)),
        ("its bottom cut along the axis: the 336 detectors (ghosted) sit in the EFCu section",
         tube(facing(ph)) + [make_actor([placed(g, i) for i in dets], GHOST, 0)[0]], (th, ph, 1.0, _bounds(pts, zmax=1500))),
    ])


def uglar_png(out):
    """Left: the underground argon alone. Right: its bottom, ghosted, with what sits in it: detectors, fibres, SiPMs."""
    g, colours = load_geometry(), read_colours()
    names = [str(n) for n in g["names"]]
    u = names.index(UGLAR)
    dets = [i for i, n in enumerate(names) if re.fullmatch(r"V\d{5}Z", n)]
    fibres = [i for i, n in enumerate(names) if n.startswith("fiber_core")]
    sipms = [i for i, n in enumerate(names) if re.fullmatch(r"larinstr_.*_sipm_silicon_.*", n)]
    mounts = [i for i, n in enumerate(names) if re.fullmatch(r"larinstr_.*_sipm_wrap_copper_.*", n)]
    p = parts_info([UGLAR]).get(UGLAR, {})
    print(f"{UGLAR}: {p.get('material', '?')}, {p.get('unit_volume', 0) / 1e6:.2f} m3, {p.get('total_mass', 0) / 1e6:.2f} t; "
          f"inside it {len(dets)} detectors, {int(sum(g['n_pl'][i] for i in fibres))} fibres, {int(sum(g['n_pl'][i] for i in sipms))} SiPMs")
    inner = np.concatenate([placed(g, i)[0] for i in dets + fibres])
    th, ph = START_THETA_PHI
    write_panels(out, [
        ("the underground argon (UGLAr), alone: the argon veto of sim/run.mac", [make_actor([placed(g, u)], colours[UGLAR], 0)[0]],
         (th, ph, 0.9)),
        ("its bottom, ghosted: the 336 detectors, and the fibres whose SiPMs see its light",
         [make_actor([placed(g, u)], (0.55, 0.8, 1.0, 0.18), 0)[0], make_actor([placed(g, i) for i in dets], colours[names[dets[0]]], 0)[0],
          make_actor([placed(g, i) for i in fibres], (*colours[names[fibres[0]]][:3], 0.22), 0)[0],  # ghosted: 12096 fibres hide all else
          make_actor([placed(g, i) for i in sipms], colours[names[sipms[0]]], 0)[0],
          make_actor([placed(g, i) for i in mounts], colours[names[mounts[0]]], 0)[0]],
         (th, ph, 1.0, _bounds(inner))),
    ])


def hpge_png(out):
    """Left: all the HPGe detectors alone. Right: one of them cut in half, its borehole open."""
    g, colours = load_geometry(), read_colours()
    names = [str(n) for n in g["names"]]
    dets = [i for i, n in enumerate(names) if re.fullmatch(r"V\d{5}Z", n)]
    d = names.index(UNIT)
    p = parts_info([names[i] for i in dets])
    mass = sum(p[names[i]]["total_mass"] for i in dets if names[i] in p)
    one = placed(g, d)[0]
    c = g["tra"][g["o_pl"][d]]
    r = np.hypot(*(one[:, :2] - c[:2]).T)
    z = one[:, 2] - c[2]
    bore = z[(r > 0.5) & (r < 6) & (z > 1)]  # the borehole wall, near the axis, above the base
    strings = len({names[i][1:4] for i in dets})
    print(f"{len(dets)} detectors in {strings} strings, {p.get(UNIT, {}).get('material', '?')}, {mass / 1e3:.1f} kg; "
          f"{UNIT}: {z.max() - z.min():.1f} mm tall, {2 * r.max():.1f} mm across, borehole {z.max() - bore.min():.0f} mm deep "
          f"from the top, {p.get(UNIT, {}).get('total_mass', 0) / 1e3:.3f} kg")
    blue = colours[UNIT]
    th, ph = START_THETA_PHI
    write_panels(out, [
        (f"all {len(dets)} HPGe detectors, alone: {strings} strings in 7 clusters of 6",
         [make_actor([placed(g, i) for i in dets], blue, 0)[0]], (28, 30, 1.1)),
        (f"one detector ({UNIT}), an ICPC, cut in half: its borehole from the top",
         [half_solid(UNIT, c, ph, blue)], (65, ph + 40, 1.0)),  # 40 deg off the cut: its face and the outside both show
    ])


PICTURES = {"--holders": holders_png, "--tube": tube_png, "--uglar": uglar_png, "--hpge": hpge_png}


def main():
    if len(sys.argv) == 3 and sys.argv[1] in PICTURES:
        return PICTURES[sys.argv[1]](Path(sys.argv[2]))
    t0 = time.monotonic()
    layers, n_tri = assemble(load_geometry(), read_colours())

    print(f"\n {'#':>2}  {'layer':<76} {'LVs':>5} {'placed':>7}")
    for j, lay in enumerate(layers):
        num = str(j + 1) if j < len(LAYERS) else "-"
        print(f" {num:>2}  {lay['label']:<76} {lay['lvs']:>5} {lay['placements']:>7}")
    n_actors = sum(len(lay["actors"]) for lay in layers)
    print(f"\n{n_tri / 1e6:.2f} M triangles in {n_actors} actors, ready in {time.monotonic() - t0:.1f} s")

    vtkMapper.SetResolveCoincidentTopologyToPolygonOffset()
    ren = vtkRenderer()
    ren.GradientBackgroundOn()
    ren.SetBackground(0.05, 0.05, 0.07)
    ren.SetBackground2(0.20, 0.22, 0.26)
    ren.SetUseDepthPeeling(True)  # for the fading layer only: opaque frames skip it
    for lay in layers:
        for actor, _ in lay["actors"]:
            ren.AddActor(actor)
    bounds = ren.ComputeVisiblePropBounds()
    start_view(ren, bounds)

    info = text_actor(15, (0.95, 0.95, 0.95), 0.012, 0.985)
    keys = text_actor(12, (0.65, 0.68, 0.72), 0.012, 0.015)
    keys.SetInput("drag rotate   wheel zoom   shift+drag pan   |   space play/pause   left/right step"
                  "   up/down speed   0 restart   c start view   r fit   f fly to cursor   q quit")
    ren.AddViewProp(info)
    ren.AddViewProp(keys)

    win = vtkRenderWindow()
    win.AddRenderer(ren)
    win.SetSize(1500, 1000)
    win.SetWindowName("LEGEND-1000 autopeel")
    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(win)
    style = vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)

    peel = Peeler(layers, win, ren, info)
    actions = {
        "space": peel.toggle,
        "Right": peel.forward,
        "Left": peel.back,
        "Up": lambda: peel.speed(1 / 1.5),
        "Down": lambda: peel.speed(1.5),
        "0": peel.restart,
        "c": lambda: start_view(ren, bounds),
        "q": iren.TerminateApp,
    }

    def on_key(_obj, _event):
        action = actions.get(iren.GetKeySym())
        if action:
            action()
            peel.dirty = True

    def on_char(_obj, _event):
        if iren.GetKeySym() in ("f", "r"):  # keep VTK's fly-to and fit; drop its other keys
            style.OnChar()

    # observers on the style replace its default handlers for these events
    style.AddObserver("KeyPressEvent", on_key)
    style.AddObserver("CharEvent", on_char)

    iren.Initialize()
    iren.AddObserver("TimerEvent", lambda *_: peel.tick())
    iren.CreateRepeatingTimer(30)
    iren.Start()


if __name__ == "__main__":
    main()
