// LEGEND.app's rt: the re-entrance tube designs, the GDML's tube in 3D. Two sliders set L1 (the SS ends at that depth
// below the tube's top) and L2 (the EFCu starts there), Cu in between; the tube's colour is that design's BI after
// all cuts over the as-built tube's, from ana/background.C's design scan (output/designs.json). It reads designs.json
// + tube.bin from its folder of the bundle (both written by make_app.py), so it needs nothing else installed.
// legend.swift opens it; --snapshot out.png --open rt [--l1 m] [--l2 m] [--size WxH] draws the window (the 3D view
// and the panel) offscreen.

import AppKit
import SceneKit
import simd

// ── the files ───────────────────────────────────────────────────────────────────────────────────────

struct Designs: Decodable {
    let title: String
    let depthTot, L1built, L2built: Double  // m below the tube's top: its length, its seams as built
    let biBuilt, tol: Double                // the as-built BI after all cuts [cts/(keV kg yr)]; the rise a design may have
    let l1, l2: [Double]                    // m: the depths the sliders step through (2 cm apart, through the seams as built)
    let mSS, mCu1: [Double]                 // kg: the SS above l1[i]; the Cu-section mass above l1[i]
    let mCu2, mEF: [Double]                 // kg: the Cu-section mass above l2[j]; the EFCu below l2[j]
    let ratio: [[Double?]]                  // the BI of design (l1[i], l2[j]) over the as-built BI; null where l1 >= l2
    let zTop: Double                        // mm: the tube's top in the GDML (z up)
    let rStep: Double                       // mm
    let radius: [Double]                    // mm: the tube's outer radius every rStep down from its top
    let nPts, nTri: Int                     // tube.bin: points and normals (float32 x3), triangles (uint32 x3)
}

func loadDesigns(_ dir: URL) throws -> (Designs, Data) {
    let d = try JSONDecoder().decode(Designs.self, from: Data(contentsOf: dir.appendingPathComponent("designs.json")))
    let bin = try Data(contentsOf: dir.appendingPathComponent("tube.bin"))
    guard bin.count == d.nPts * 24 + d.nTri * 12 else {
        throw NSError(domain: "tube.bin", code: 1, userInfo: [NSLocalizedDescriptionKey: "tube.bin does not match designs.json"])
    }
    return (d, bin)
}

func tubeGeometry(_ d: Designs, _ bin: Data) -> SCNGeometry {
    let fv = d.nPts * 12
    let source = { (offset: Int, s: SCNGeometrySource.Semantic) in
        SCNGeometrySource(data: bin.subdata(in: offset..<offset + fv), semantic: s, vectorCount: d.nPts, usesFloatComponents: true,
                          componentsPerVector: 3, bytesPerComponent: 4, dataOffset: 0, dataStride: 12)
    }
    let element = SCNGeometryElement(data: bin.subdata(in: 2 * fv..<bin.count), primitiveType: .triangles,
                                     primitiveCount: d.nTri, bytesPerIndex: 4)
    return SCNGeometry(sources: [source(0, .vertex), source(fv, .normal)], elements: [element])
}

// ── colours ───────────────────────────────────────────────────────────────────────────────────

// blue at 0.9 x the as-built BI, through green, to yellow at 1.3 and above
let stops: [(x: Double, rgb: [Double])] = [(0.9, [38, 70, 150]), (1.0, [33, 120, 190]), (1.1, [40, 170, 150]),
                                           (1.2, [140, 200, 80]), (1.3, [250, 220, 60])]

func colour(_ x: Double) -> NSColor {
    var c = stops.last!.rgb
    if x <= stops[0].x { c = stops[0].rgb }
    else if let k = stops.indices.dropFirst().first(where: { x <= stops[$0].x }) {
        let f = (x - stops[k - 1].x) / (stops[k].x - stops[k - 1].x)
        c = zip(stops[k - 1].rgb, stops[k].rgb).map { $0 + f * ($1 - $0) }
    }
    return NSColor(srgbRed: c[0] / 255, green: c[1] / 255, blue: c[2] / 255, alpha: 1)
}

let ink = NSColor(srgbRed: 0.11, green: 0.14, blue: 0.19, alpha: 1)   // text on a light tube
let passColour = NSColor(srgbRed: 0.37, green: 0.78, blue: 0.48, alpha: 1)
let failColour = NSColor(srgbRed: 1.00, green: 0.48, blue: 0.45, alpha: 1)

// the tube's inner wall, a shade darker, so its open top reads as open
let inside = """
#pragma body
float3 p = (scn_node.inverseModelViewTransform * float4(_surface.position, 1.0)).xyz;
float3 n = (scn_node.inverseModelViewTransform * float4(_surface.normal, 0.0)).xyz;
if (dot(n.xy, p.xy) < 0.0) { _surface.diffuse.rgb *= 0.4; }
"""

func flat(_ c: NSColor, alpha: CGFloat = 1) -> SCNMaterial {
    let m = SCNMaterial()
    m.lightingModel = .constant
    m.diffuse.contents = c
    m.transparency = alpha
    m.isDoubleSided = true
    return m
}

// ── the scene ───────────────────────────────────────────────────────────────────────────────────────

func label3D(_ size: CGFloat) -> SCNNode {
    let t = SCNText(string: "", extrusionDepth: 0)
    t.font = .systemFont(ofSize: size, weight: .semibold)
    t.firstMaterial = flat(fg)
    return SCNNode(geometry: t)
}

// the text, its anchor (align 0 left, 0.5 centre, 1 right; middle of the line) at p
func place(_ node: SCNNode, _ s: String, _ p: SCNVector3, align: CGFloat, colour: NSColor) {
    let t = node.geometry as! SCNText
    if (t.string as? String) != s { t.string = s }
    t.firstMaterial!.diffuse.contents = colour
    node.isHidden = s.isEmpty
    let (a, b) = node.boundingBox
    node.pivot = SCNMatrix4MakeTranslation(a.x + (b.x - a.x) * align, (a.y + b.y) / 2, 0)
    node.position = p
}

final class TubeScene {
    let d: Designs
    let scene = SCNScene()
    let camera = SCNNode()
    let body = SCNNode()                            // the tube
    let paint = SCNMaterial()                       // its colour: the design's BI
    let seams = [SCNNode(), SCNNode()]              // L1 and L2: white rings
    let rig = SCNNode()                             // the labels: it turns about the vertical to face the camera
    let sections = [label3D(150), label3D(150), label3D(150)]   // SS, Cu, EFCu, on the tube's front
    let seamLabels = [label3D(110), label3D(110)]   // right of the rings
    var target = SCNVector3(0, 0, 0)

    // the scene's y is the GDML's z; depths are m below the tube's top
    func y(_ depth: Double) -> CGFloat { CGFloat(d.zTop - depth * 1000) }
    func r(_ depth: Double) -> CGFloat {
        CGFloat(d.radius[min(max(Int((depth * 1000 / d.rStep).rounded()), 0), d.radius.count - 1)])
    }

    init(_ d: Designs, _ bin: Data) {
        self.d = d
        let world = SCNNode()
        world.eulerAngles.x = -.pi / 2   // the GDML is z up, SceneKit y up
        let geometry = tubeGeometry(d, bin)
        paint.lightingModel = .lambert
        paint.isDoubleSided = true
        paint.shaderModifiers = [.surface: inside]
        geometry.materials = [paint]
        body.geometry = geometry
        world.addChildNode(body)
        scene.rootNode.addChildNode(world)

        for depth in [d.L1built, d.L2built] { scene.rootNode.addChildNode(dashedRing(depth)) }
        for s in seams {
            let t = SCNTorus(ringRadius: 1000, pipeRadius: 12)
            t.ringSegmentCount = 144
            t.pipeSegmentCount = 12
            t.firstMaterial = flat(fg)
            s.geometry = t
            scene.rootNode.addChildNode(s)
        }
        let face = SCNBillboardConstraint()
        face.freeAxes = .Y
        rig.constraints = [face]
        for depth in [d.L1built, d.L2built] {
            let t = label3D(80)
            place(t, "as built", SCNVector3(-(r(depth) + 60), y(depth), 0), align: 1, colour: muted)
            rig.addChildNode(t)
        }
        for t in sections + seamLabels { rig.addChildNode(t) }
        scene.rootNode.addChildNode(rig)

        scene.background.contents = gradient(bottom: NSColor(srgbRed: 0.05, green: 0.05, blue: 0.07, alpha: 1),
                                             top: NSColor(srgbRed: 0.20, green: 0.22, blue: 0.26, alpha: 1))
        let ambient = SCNNode()
        ambient.light = SCNLight()
        ambient.light!.type = .ambient
        ambient.light!.intensity = 110
        scene.rootNode.addChildNode(ambient)
        camera.camera = SCNCamera()
        camera.camera!.fieldOfView = 30
        camera.camera!.automaticallyAdjustsZRange = true
        let head = SCNNode(), key = SCNNode()   // a headlight, and a key light from the upper left that rounds the tube
        head.light = SCNLight()
        head.light!.type = .directional
        head.light!.intensity = 330
        key.light = SCNLight()
        key.light!.type = .directional
        key.light!.intensity = 700
        key.eulerAngles = SCNVector3(-0.45, -0.75, 0)
        camera.addChildNode(head)
        camera.addChildNode(key)
        scene.rootNode.addChildNode(camera)
        startView()
    }

    // the seam as built: a dashed ring, faint
    func dashedRing(_ depth: Double) -> SCNNode {
        let ring = SCNNode(), n = 48, rr = r(depth) + 6
        let dash = SCNBox(width: 2 * .pi * rr / CGFloat(n) / 2, height: 10, length: 10, chamferRadius: 0)
        dash.firstMaterial = flat(fg, alpha: 0.55)
        for k in 0..<n {
            let turn = SCNNode(), piece = SCNNode(geometry: dash)
            turn.eulerAngles.y = 2 * .pi * CGFloat(k) / CGFloat(n)
            piece.position = SCNVector3(0, 0, rr)   // its length along the ring
            turn.addChildNode(piece)
            ring.addChildNode(turn)
        }
        ring.position = SCNVector3(0, y(depth), 0)
        return ring
    }

    // the design: SS down to L1, Cu down to L2, EFCu below; the tube in the colour of its BI ratio
    func show(L1: Double, L2: Double, ratio: Double) {
        let c = colour(ratio)
        paint.diffuse.contents = c
        for (s, depth) in zip(seams, [L1, L2]) {
            (s.geometry as! SCNTorus).ringRadius = r(depth) + 4
            s.position = SCNVector3(0, y(depth), 0)
        }
        let rgb = c.usingColorSpace(.sRGB)!
        let light = 0.2126 * rgb.redComponent + 0.7152 * rgb.greenComponent + 0.0722 * rgb.blueComponent > 0.5
        let edges = [0, L1, L2, d.depthTot], names = ["SS", "Cu", "EFCu"]
        for k in 0..<3 {
            let mid = (edges[k] + edges[k + 1]) / 2
            place(sections[k], edges[k + 1] - edges[k] > 0.25 ? names[k] : "", SCNVector3(0, y(mid), r(mid) + 40),
                  align: 0.5, colour: light ? ink : fg)
        }
        for (k, depth) in [L1, L2].enumerated() {
            place(seamLabels[k], String(format: "L%d  %.2f m", k + 1, depth), SCNVector3(r(depth) + 60, y(depth), 0),
                  align: 0, colour: fg)
        }
    }

    // the start view: from 10 deg above the horizontal, framing the tube
    func startView() {
        let th = 80.0 * .pi / 180, ph = 20.0 * .pi / 180
        frame(direction: SCNVector3(sin(th) * cos(ph), cos(th), -sin(th) * sin(ph)))   // (x, y, z) z up -> (x, z, -y)
    }

    func frame(direction dir: SCNVector3) {
        let (a, b) = body.boundingBox
        var lo = SIMD3<Double>(repeating: .infinity), hi = SIMD3<Double>(repeating: -.infinity)
        for c in [a, b] {
            let w = body.convertPosition(c, to: nil)
            lo = simd_min(lo, SIMD3(Double(w.x), Double(w.y), Double(w.z)))
            hi = simd_max(hi, SIMD3(Double(w.x), Double(w.y), Double(w.z)))
        }
        let centre = (lo + hi) / 2, radius = (hi.y - lo.y) / 2 * 1.12   // the tube's height, with room for the labels
        let d = simd_normalize(SIMD3(Double(dir.x), Double(dir.y), Double(dir.z)))
        let eye = centre + d * radius / sin(Double(camera.camera!.fieldOfView) * .pi / 360)
        camera.position = SCNVector3(eye.x, eye.y, eye.z)
        camera.look(at: SCNVector3(centre.x, centre.y, centre.z), up: SCNVector3(0, 1, 0), localFront: SCNVector3(0, 0, -1))
        target = SCNVector3(centre.x, centre.y, centre.z)
    }
}

// ── the panel: the colour scale ────────────────────────────────────────────────────────────────

final class Scale: NSView {   // the colour scale, 0.9 to 1.3
    override func draw(_ dirty: NSRect) {
        let at = stops.map { CGFloat(($0.x - stops[0].x) / (stops.last!.x - stops[0].x)) }
        NSGradient(colors: stops.map { colour($0.x) }, atLocations: at, colorSpace: .sRGB)!
            .draw(in: NSBezierPath(roundedRect: bounds, xRadius: bounds.height / 2, yRadius: bounds.height / 2), angle: 0)
    }
}

// ── the window ──────────────────────────────────────────────────────────────────────────────────────

func nearest(_ a: [Double], _ v: Double) -> Int { a.indices.min { abs(a[$0] - v) < abs(a[$1] - v) }! }

final class RTController: NSObject, Viewer {
    let d: Designs
    let tube: TubeScene
    let window: NSWindow
    let view: KeyView
    let panel: Fill
    var sliders: [NSSlider] = []
    let depthText = [field(13, .regular, muted), field(13, .regular, muted)]
    let ratioText = field(34, .bold)
    let verdictText = field(13, .semibold)
    let massText = [field(13), field(13), field(13)]
    let built: (Int, Int)   // the slider steps nearest the seams as built
    static let P: CGFloat = 360   // the panel's width
    let noteText = wrapping("", 12, muted, width: RTController.P - 72)

    init(_ dir: URL) throws {
        let (d, bin) = try loadDesigns(dir)
        let t0 = Date()
        self.d = d
        tube = TubeScene(d, bin)
        built = (nearest(d.l1, d.L1built), nearest(d.l2, d.L2built))
        let screen = NSScreen.main?.visibleFrame ?? NSRect(x: 0, y: 0, width: 1500, height: 1000)
        let size = NSSize(width: min(1500, screen.width * 0.9), height: min(1000, screen.height * 0.9)), P = Self.P
        window = NSWindow(contentRect: NSRect(origin: .zero, size: size), styleMask: [.titled, .closable, .miniaturizable, .resizable],
                          backing: .buffered, defer: false)
        window.isReleasedWhenClosed = false
        window.title = d.title
        window.appearance = NSAppearance(named: .darkAqua)
        window.minSize = NSSize(width: 800, height: 720)
        view = KeyView(frame: NSRect(x: 0, y: 0, width: size.width - P, height: size.height))
        panel = Fill(frame: NSRect(x: size.width - P, y: 0, width: P, height: size.height))
        super.init()

        view.autoresizingMask = [.width, .height]
        view.scene = tube.scene
        view.pointOfView = tube.camera
        view.allowsCameraControl = true
        view.defaultCameraController.interactionMode = .orbitTurntable
        view.defaultCameraController.worldUp = SCNVector3(0, 1, 0)
        view.defaultCameraController.target = tube.target
        view.antialiasingMode = .multisampling4X
        window.contentView!.addSubview(view)
        let (tube, view) = (self.tube, self.view)
        view.keys = ["c": { tube.startView(); view.defaultCameraController.target = tube.target },
                     "0": { [unowned self] in self.asBuilt() }, "q": { NSApp.terminate(nil) }]

        panel.autoresizingMask = [.minXMargin, .height]
        window.contentView!.addSubview(panel)
        let W = P - 40
        sliders = [d.l1.count, d.l2.count].map {
            let s = NSSlider(value: 0, minValue: 0, maxValue: Double($0 - 1), target: self, action: #selector(slid(_:)))
            s.isContinuous = true
            return s
        }
        let reset = NSButton(title: "Back to as built", target: self, action: #selector(asBuilt))
        reset.bezelStyle = .rounded
        let ticks = NSStackView(views: ["0.9", "1.0", "1.1", "1.2", "≥ 1.3"].map { field(11, .regular, muted, $0) })
        ticks.distribution = .equalSpacing
        let scale = Scale()
        scale.heightAnchor.constraint(equalToConstant: 12).isActive = true
        let stack = column([
            field(19, .bold, fg, "Re-entrance tube designs"),
            wrapping("SS down to L1, Cu down to L2, EFCu below. The tube's colour is that design's BI after all cuts "
                     + "over the as-built tube's.", 13, muted, width: W),
            card([row(field(13, .semibold, fg, "L1: SS ends at depth"), depthText[0]), sliders[0],
                  row(field(13, .semibold, fg, "L2: EFCu starts at depth"), depthText[1]), sliders[1], row(reset, NSView())], width: W),
            card([field(13, .regular, muted, "BI over the as-built tube"), ratioText, scale, ticks,
                  row(field(13, .regular, fg, "within 10% of the as-built BI?"), verdictText)], width: W),
            card([row(field(13, .semibold, fg, "SS"), massText[0]), row(field(13, .semibold, fg, "Cu"), massText[1]),
                  row(field(13, .semibold, fg, "EFCu"), massText[2]), noteText], width: W),
        ], spacing: 14, width: W)
        let keys = wrapping("drag rotate   scroll zoom   c start view   0 as built   q quit", 11, muted, width: W)
        keys.translatesAutoresizingMaskIntoConstraints = false
        panel.addSubview(stack)
        panel.addSubview(keys)
        NSLayoutConstraint.activate([
            stack.leadingAnchor.constraint(equalTo: panel.leadingAnchor, constant: 20),
            stack.topAnchor.constraint(equalTo: panel.topAnchor, constant: 20),
            keys.leadingAnchor.constraint(equalTo: panel.leadingAnchor, constant: 20),
            keys.bottomAnchor.constraint(equalTo: panel.bottomAnchor, constant: -14),
        ])
        select(built.0, built.1)
        print(String(format: "%@: %d designs, the tube in %.2f M triangles, ready in %.1f s", d.title,
                     d.ratio.joined().compactMap { $0 }.count, Double(d.nTri) / 1e6, Date().timeIntervalSince(t0)))
    }

    func show() {
        window.center()
        window.makeKeyAndOrderFront(nil)
        window.makeFirstResponder(view)
    }

    func close() { view.keys = [:] }

    @objc func slid(_ s: NSSlider) { update(moved: s === sliders[0] ? 0 : 1) }

    @objc func asBuilt() { select(built.0, built.1) }

    func select(_ i: Int, _ j: Int) {
        sliders[0].doubleValue = Double(i)
        sliders[1].doubleValue = Double(j)
        update(moved: 1)
    }

    // keep some Cu between the seams, then show the design
    func update(moved: Int) {
        var i = Int(sliders[0].doubleValue.rounded()), j = Int(sliders[1].doubleValue.rounded())
        if d.l1[i] >= d.l2[j] {
            if moved == 0 { j = d.l2.firstIndex { $0 > d.l1[i] } ?? d.l2.count - 1 }
            if d.l1[i] >= d.l2[j] { i = d.l1.lastIndex { $0 < d.l2[j] } ?? 0 }
        }
        sliders[0].doubleValue = Double(i)
        sliders[1].doubleValue = Double(j)
        let L1 = d.l1[i], L2 = d.l2[j], c = d.ratio[i][j] ?? .nan
        tube.show(L1: L1, L2: L2, ratio: c)
        depthText[0].stringValue = String(format: "%.2f m", L1)
        depthText[1].stringValue = String(format: "%.2f m", L2)
        ratioText.stringValue = String(format: c < 10 ? "× %.2f" : "× %.0f", c)
        let ok = c <= 1 + d.tol + 1e-9
        verdictText.stringValue = String(format: "%@%.0f%% · %@", c >= 1 ? "+" : "", 100 * (c - 1), ok ? "passes" : "fails")
        verdictText.textColor = ok ? passColour : failColour
        let ss = d.mSS[i], cu = d.mCu2[j] - d.mCu1[i], ef = d.mEF[j], tot = ss + cu + ef
        for (t, m) in zip(massText, [ss, cu, ef]) { t.stringValue = String(format: "%.0f kg · %.0f%%", m, 100 * m / tot) }
        noteText.stringValue = (i, j) == built
            ? String(format: "As built. BI after all cuts %.2e cts/(keV·kg·yr).", d.biBuilt)
            : String(format: "A design passes when its BI rises by at most %.0f%%. Dashed rings: the seams as built.", 100 * d.tol)
    }

    // --snapshot out.png [--l1 m] [--l2 m] [--size WxH]: the 3D view and the panel side by side, as the window shows them
    func snapshot(to path: String, _ args: [String]) {
        let size = argument(args, "--size").flatMap(parseSize) ?? NSSize(width: 1500, height: 1000), P = Self.P
        window.setContentSize(size)
        select(argument(args, "--l1").flatMap(Double.init).map { nearest(d.l1, $0) } ?? built.0,
               argument(args, "--l2").flatMap(Double.init).map { nearest(d.l2, $0) } ?? built.1)
        window.contentView!.layoutSubtreeIfNeeded()
        guard let rep = panel.bitmapImageRepForCachingDisplay(in: panel.bounds) else { return }
        panel.cacheDisplay(in: panel.bounds, to: rep)
        let r = SCNRenderer(device: MTLCreateSystemDefaultDevice(), options: nil)
        r.scene = tube.scene
        r.pointOfView = tube.camera
        let w = size.width - P
        let scene = r.snapshot(atTime: 0, with: CGSize(width: w, height: size.height), antialiasingMode: .multisampling4X)
        let image = NSImage(size: size)
        image.lockFocus()
        scene.draw(in: NSRect(x: 0, y: 0, width: w, height: size.height))
        rep.draw(in: NSRect(x: w, y: 0, width: P, height: size.height))
        image.unlockFocus()
        writePNG(image, path)
        print("wrote \(path): L1 \(depthText[0].stringValue), L2 \(depthText[1].stringValue), \(ratioText.stringValue)")
    }
}
