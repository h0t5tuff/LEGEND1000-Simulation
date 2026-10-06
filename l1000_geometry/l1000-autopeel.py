#!/usr/bin/env python
"""Peel the LEGEND-1000 geometry layer by layer on a timer while you rotate it (see README.md).

    ~/venvs/v/bin/python l1000-autopeel.py
"""

import itertools
import re
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
CACHE = HERE / "output" / "autopeel-cache.npz"
DWELL = 3.0  # s per layer at start (up/down change it)
FADE = 0.5  # s a layer takes to fade out or back in
START_THETA_PHI = (70, 20)  # deg, z up: the start view of l1000-vis.mac

# Peel order, outside in. Each LV goes to the first layer whose regex matches its whole name;
# whatever no layer matches (the HPGe detectors V*) is left at the end.
LAYERS = [
    ("outer cryostat: outer steel tank, skirt and foot", r"outercryostat|skirt|foot"),
    ("vacuum gap between the two steel tanks", r"vacuumgap"),
    ("inner cryostat", r"innercryostat"),
    ("atmospheric LAr (outer bath)", r"atmosphericlar"),
    ("neutron moderator", r"neutronmoderator"),
    ("outer reflector: TPB coating", r"wls_tpb_outer_atmospheric_lv"),
    ("outer reflector: tetratex lining", r"wls_tetratex_outer_atmospheric_lv"),
    ("re-entrant tube: EFCu (z < 852 mm), OFHC Cu (852-2857), SS 316L (> 2857)",
     r"reentrancetube|ofhc_cu|ss_316l"),
    ("underground LAr (inner active bath)", r"undergroundlar"),
    ("inner reflector: TPB coating", r"wls_tpb_inner_argon_lv"),
    ("inner reflector: tetratex lining", r"wls_tetratex_inner_argon_lv"),
    ("WLS fibres: outer and inner cladding, TPB coating, PS core", r"fiber_.*"),
    ("SiPM mounts (copper rings, top and bottom)", r"sipm_outer_.*"),
    ("SiPMs", r"sipm_r.*"),
    ("string hardware: rods, weldments, support structure, signal and HV cables",
     r"hpge_support_copper_(rod_string_\d+|string_support_structure|weldment_top_V\d+)"
     r"|cable_(hv|signal)_V\d+"),
    ("ultem insulators and signal/HV clamps", r"ultem_.*"),
    ("PEN holder plates", r"pen_xlarge|hpge_support_copper_tristar_xlarge"),
    ("signal ASICs (1 mm silica chips under each detector: look from below)", r"signal_asic_V\d+"),
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


def placed(g, i):
    """LV i's surface at all of its placements, in world coordinates."""
    p = g["pts"][g["o_pts"][i]:g["o_pts"][i + 1]]
    n = g["nrm"][g["o_pts"][i]:g["o_pts"][i + 1]]
    t = g["tri"][g["o_tri"][i]:g["o_tri"][i + 1]]
    rot = g["rot"][g["o_pl"][i]:g["o_pl"][i + 1]]
    tra = g["tra"][g["o_pl"][i]:g["o_pl"][i + 1]]
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


def make_actor(parts, rgba, depth):
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
    mapper.SetInputData(pd)
    mapper.ScalarVisibilityOff()
    mapper.StaticOn()  # geometry never changes: upload once
    # a daughter that shares a face with its mother (e.g. ofhc_cu / ss_316l on the re-entrant
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

    core = (f"the {len(left)} HPGe detectors" if all(re.fullmatch(r"V\d+", n) for n in left)
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


def main():
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
