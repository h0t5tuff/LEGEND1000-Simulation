// The autopeel as a native macOS app: peels a geometry layer by layer on a timer while you rotate it, like
// LEGEND200/l200-autopeel.py and geom/l1000-autopeel.py. It reads geometry.json + geometry.bin from its own
// bundle (both written by make_app.py), so it needs nothing else installed. Built by make_app.py:
//   xcrun swiftc -O -parse-as-library -swift-version 5 autopeel.swift
// Hidden check: <app>/Contents/MacOS/<name> --snapshot out.png [--peel N] renders one frame offscreen and quits.

import AppKit
import SceneKit
import simd

// ── the geometry file ───────────────────────────────────────────────────────────────────────────────

struct Part: Decodable { let name: String; let pts: [Int]; let tri: [Int]; let pl: [Int] }  // [first, count]
struct Group: Decodable { let rgba: [Float]; let depth: Int; let parts: [Int] }
struct Layer: Decodable { let label: String; let groups: [Group] }
struct Geometry: Decodable {
    let title: String
    let start_theta_phi: [Double]   // deg, z up
    let dwell: Double                // s per layer at start
    let fade: Double                 // s a layer takes to fade
    let layers: [Layer]              // outside in; the last one is the core, which stays
    let parts: [Part]
    let arrays: [String: [Int]]      // pts, nrm (float32 x3), tri (uint32 x3), rot (float32 x9), tra (float32 x3): [byte offset, count]
}

func loadGeometry(_ dir: URL) throws -> (Geometry, Data) {
    let geo = try JSONDecoder().decode(Geometry.self, from: Data(contentsOf: dir.appendingPathComponent("geometry.json")))
    return (geo, try Data(contentsOf: dir.appendingPathComponent("geometry.bin"), options: .alwaysMapped))
}

func array<T>(_ bin: Data, _ geo: Geometry, _ key: String, _ type: T.Type) -> [T] {
    let (offset, count) = (geo.arrays[key]![0], geo.arrays[key]![1])
    return bin.withUnsafeBytes { raw in
        Array(UnsafeRawBufferPointer(rebasing: raw[offset..<offset + count * MemoryLayout<T>.stride]).bindMemory(to: T.self))
    }
}

// one mesh per (layer, colour, depth) group: every part of the group at all of its placements, in world coordinates
func groupGeometry(_ g: Group, _ geo: Geometry, pts: [Float], nrm: [Float], tri: [UInt32], rot: [Float], tra: [Float]) -> SCNGeometry {
    var nv = 0, nt = 0
    for i in g.parts { nv += geo.parts[i].pts[1] * geo.parts[i].pl[1]; nt += geo.parts[i].tri[1] * geo.parts[i].pl[1] }
    var pos = [Float](repeating: 0, count: nv * 3), nor = [Float](repeating: 0, count: nv * 3)
    var idx = [UInt32](repeating: 0, count: nt * 3)
    var v = 0, t = 0
    for i in g.parts {
        let p = geo.parts[i]
        for k in p.pl[0]..<(p.pl[0] + p.pl[1]) {
            let r = rot[(k * 9)..<(k * 9 + 9)].map { $0 }
            let R = simd_float3x3(rows: [SIMD3(r[0], r[1], r[2]), SIMD3(r[3], r[4], r[5]), SIMD3(r[6], r[7], r[8])])
            let N = R.inverse.transpose
            let T = SIMD3(tra[k * 3], tra[k * 3 + 1], tra[k * 3 + 2])
            for a in 0..<p.pts[1] {
                let j = (p.pts[0] + a) * 3, o = (v + a) * 3
                let q = R * SIMD3(pts[j], pts[j + 1], pts[j + 2]) + T
                let n = simd_normalize(N * SIMD3(nrm[j], nrm[j + 1], nrm[j + 2]))
                pos[o] = q.x; pos[o + 1] = q.y; pos[o + 2] = q.z
                nor[o] = n.x; nor[o + 1] = n.y; nor[o + 2] = n.z
            }
            for b in 0..<(p.tri[1] * 3) { idx[t * 3 + b] = UInt32(v) + tri[p.tri[0] * 3 + b] }
            v += p.pts[1]; t += p.tri[1]
        }
    }
    let source = { (a: [Float], s: SCNGeometrySource.Semantic) in
        SCNGeometrySource(data: a.withUnsafeBufferPointer { Data(buffer: $0) }, semantic: s, vectorCount: nv,
                          usesFloatComponents: true, componentsPerVector: 3, bytesPerComponent: 4, dataOffset: 0, dataStride: 12)
    }
    let element = SCNGeometryElement(data: idx.withUnsafeBufferPointer { Data(buffer: $0) }, primitiveType: .triangles,
                                     primitiveCount: nt, bytesPerIndex: 4)
    return SCNGeometry(sources: [source(pos, .vertex), source(nor, .normal)], elements: [element])
}

// a daughter sharing a face with its mother wins the depth test: its vertices move a hair towards the eye
let depthPull = """
#pragma arguments
float depthPull;
#pragma body
float4 v = scn_node.modelViewTransform * _geometry.position;
v.xyz *= (1.0 - depthPull);
_geometry.position = scn_node.inverseModelViewTransform * v;
"""

func material(_ rgba: [Float], depth: Int) -> SCNMaterial {
    let m = SCNMaterial()
    m.lightingModel = .blinn
    m.diffuse.contents = NSColor(srgbRed: CGFloat(rgba[0]), green: CGFloat(rgba[1]), blue: CGFloat(rgba[2]), alpha: 1)
    m.specular.contents = NSColor(white: 0.18, alpha: 1)
    m.shininess = 25
    m.transparency = CGFloat(rgba[3])
    m.isDoubleSided = true
    m.shaderModifiers = [.geometry: depthPull]
    m.setValue(NSNumber(value: Float(depth) * 2e-5), forKey: "depthPull")
    return m
}

// ── the scene ───────────────────────────────────────────────────────────────────────────────────────

final class Peel {
    let geo: Geometry
    let scene = SCNScene()
    let camera = SCNNode()
    var layerNodes: [SCNNode] = []   // one per layer, the core last
    var triangles = 0

    init(_ geo: Geometry, _ bin: Data) {
        self.geo = geo
        let pts = array(bin, geo, "pts", Float.self), nrm = array(bin, geo, "nrm", Float.self)
        let tri = array(bin, geo, "tri", UInt32.self), rot = array(bin, geo, "rot", Float.self), tra = array(bin, geo, "tra", Float.self)
        let world = SCNNode()
        world.eulerAngles.x = -.pi / 2   // the geometry is z up, SceneKit y up
        scene.rootNode.addChildNode(world)
        for layer in geo.layers {
            let node = SCNNode()
            for g in layer.groups {
                let geometry = groupGeometry(g, geo, pts: pts, nrm: nrm, tri: tri, rot: rot, tra: tra)
                geometry.materials = [material(g.rgba, depth: g.depth)]
                triangles += geometry.elements.reduce(0) { $0 + $1.primitiveCount }
                node.addChildNode(SCNNode(geometry: geometry))
            }
            world.addChildNode(node)
            layerNodes.append(node)
        }
        scene.background.contents = gradient(bottom: NSColor(srgbRed: 0.05, green: 0.05, blue: 0.07, alpha: 1),
                                             top: NSColor(srgbRed: 0.20, green: 0.22, blue: 0.26, alpha: 1))
        let ambient = SCNNode()
        ambient.light = SCNLight()
        ambient.light!.type = .ambient
        ambient.light!.intensity = 260
        scene.rootNode.addChildNode(ambient)
        camera.camera = SCNCamera()
        camera.camera!.fieldOfView = 30
        camera.camera!.automaticallyAdjustsZRange = true
        let head = SCNNode()   // a headlight, as in VTK
        head.light = SCNLight()
        head.light!.type = .directional
        head.light!.intensity = 850
        camera.addChildNode(head)
        scene.rootNode.addChildNode(camera)
        startView()
    }

    // the start view: theta/phi from z up, framing every layer
    func startView() {
        let th = geo.start_theta_phi[0] * .pi / 180, ph = geo.start_theta_phi[1] * .pi / 180
        let dir = SCNVector3(sin(th) * cos(ph), cos(th), -sin(th) * sin(ph))   // (x, y, z) z up -> (x, z, -y)
        frame(nodes: layerNodes, direction: dir)
    }

    // keep the direction of view, frame the given (visible) layers
    func frame(nodes: [SCNNode], direction: SCNVector3? = nil) {
        var lo = SIMD3<Double>(repeating: .infinity), hi = SIMD3<Double>(repeating: -.infinity)
        for n in nodes where !n.isHidden {
            let (a, b) = n.boundingBox
            if a.x > b.x { continue }
            for c in [a, b] {
                let w = n.convertPosition(c, to: nil)
                lo = simd_min(lo, SIMD3(Double(w.x), Double(w.y), Double(w.z)))
                hi = simd_max(hi, SIMD3(Double(w.x), Double(w.y), Double(w.z)))
            }
        }
        if lo.x > hi.x { return }
        let centre = (lo + hi) / 2, radius = simd_length(hi - lo) / 2
        var d = SIMD3<Double>(0, 0, 1)
        if let dir = direction { d = simd_normalize(SIMD3(Double(dir.x), Double(dir.y), Double(dir.z))) }
        else { let f = camera.worldFront; d = simd_normalize(-SIMD3(Double(f.x), Double(f.y), Double(f.z))) }
        let dist = radius / sin(Double(camera.camera!.fieldOfView) * .pi / 360) * 1.02
        let eye = centre + d * dist
        camera.position = SCNVector3(eye.x, eye.y, eye.z)
        camera.look(at: SCNVector3(centre.x, centre.y, centre.z), up: SCNVector3(0, 1, 0), localFront: SCNVector3(0, 0, -1))
        target = SCNVector3(centre.x, centre.y, centre.z)
    }
    var target = SCNVector3(0, 0, 0)
}

func gradient(bottom: NSColor, top: NSColor) -> NSImage {
    let image = NSImage(size: NSSize(width: 8, height: 512))
    image.lockFocus()
    NSGradient(starting: bottom, ending: top)!.draw(in: NSRect(x: 0, y: 0, width: 8, height: 512), angle: 90)
    image.unlockFocus()
    return image
}

// ── the peel: a timer that only changes visibility and opacity ─────────────────────────────────────────

final class Peeler {
    let peel: Peel
    let n: Int
    var peeled = 0, playing = true
    var dwell: Double
    var nextAt: Date
    var onChange: () -> Void = {}

    init(_ peel: Peel) {
        self.peel = peel
        n = peel.layerNodes.count - 1
        dwell = peel.geo.dwell
        nextAt = Date().addingTimeInterval(dwell)
    }

    func tick() { if playing && Date() >= nextAt { forward() } }

    func fade(_ i: Int, out: Bool) {
        let node = peel.layerNodes[i]
        node.removeAction(forKey: "fade")
        if !out { node.isHidden = false }
        let action = SCNAction.fadeOpacity(to: out ? 0 : 1, duration: peel.geo.fade * Double(out ? node.opacity : 1 - node.opacity))
        node.runAction(out ? .sequence([action, .hide()]) : action, forKey: "fade")
    }

    func forward() {
        if peeled < n { fade(peeled, out: true); peeled += 1 }
        playing = playing && peeled < n
        nextAt = Date().addingTimeInterval(dwell)
        onChange()
    }

    func back() {
        if peeled > 0 { peeled -= 1; fade(peeled, out: false) }
        nextAt = Date().addingTimeInterval(dwell)
        onChange()
    }

    func restart() {
        for i in 0..<peeled { fade(i, out: false) }
        peeled = 0; playing = true
        nextAt = Date().addingTimeInterval(dwell)
        onChange()
    }

    func toggle() {
        if peeled == n { return restart() }
        playing.toggle()
        nextAt = Date().addingTimeInterval(dwell)
        onChange()
    }

    func speed(_ factor: Double) {
        dwell = min(max(dwell * factor, 0.3), 60)
        nextAt = min(nextAt, Date().addingTimeInterval(dwell))
        onChange()
    }

    var text: String {
        var state = playing ? String(format: "peeling every %.1f s", dwell) : "paused"
        if peeled == n { state = "done (space starts again)" }
        var lines = ["\(peeled)/\(n) layers removed    \(state)"]
        if peeled > 0 { lines.append("removed:  \(peel.geo.layers[peeled - 1].label)") }
        lines.append("\(peeled < n ? "next" : "left"):  \(peel.geo.layers[peeled].label)")
        return lines.joined(separator: "\n")
    }
}

// ── the window ──────────────────────────────────────────────────────────────────────────────────────

final class PeelView: SCNView {
    var keys: [String: () -> Void] = [:]
    override var acceptsFirstResponder: Bool { true }
    override func keyDown(with event: NSEvent) {
        let name: String
        switch event.keyCode {
        case 49: name = "space"
        case 123: name = "Left"
        case 124: name = "Right"
        case 125: name = "Down"
        case 126: name = "Up"
        default: name = event.charactersIgnoringModifiers ?? ""
        }
        if let f = keys[name] { f() } else { super.keyDown(with: event) }
    }
}

func label(_ size: CGFloat, _ colour: NSColor) -> NSTextField {
    let t = NSTextField(labelWithString: "")
    t.font = .systemFont(ofSize: size)
    t.textColor = colour
    t.maximumNumberOfLines = 0
    let shadow = NSShadow()
    shadow.shadowColor = .black
    shadow.shadowOffset = NSSize(width: 1, height: -1)
    shadow.shadowBlurRadius = 2
    t.shadow = shadow
    t.translatesAutoresizingMaskIntoConstraints = false
    return t
}

final class AppDelegate: NSObject, NSApplicationDelegate {
    var window: NSWindow!
    var peel: Peel!
    var peeler: Peeler!
    var timer: Timer?

    func applicationDidFinishLaunching(_ note: Notification) {
        let (geo, bin): (Geometry, Data)
        do { (geo, bin) = try loadGeometry(Bundle.main.resourceURL!) } catch {
            let a = NSAlert(); a.messageText = "No geometry in the app"; a.informativeText = "\(error)"; a.runModal()
            NSApp.terminate(nil); return
        }
        let t0 = Date()
        peel = Peel(geo, bin)
        peeler = Peeler(peel)
        let args = CommandLine.arguments
        if let i = args.firstIndex(of: "--snapshot"), i + 1 < args.count { return snapshot(to: args[i + 1], args) }

        let screen = NSScreen.main?.visibleFrame ?? NSRect(x: 0, y: 0, width: 1500, height: 1000)
        let size = NSSize(width: min(1500, screen.width * 0.9), height: min(1000, screen.height * 0.9))
        window = NSWindow(contentRect: NSRect(origin: .zero, size: size), styleMask: [.titled, .closable, .miniaturizable, .resizable],
                          backing: .buffered, defer: false)
        window.title = "\(geo.title) autopeel"
        window.center()
        let view = PeelView(frame: window.contentView!.bounds)
        view.autoresizingMask = [.width, .height]
        view.scene = peel.scene
        view.pointOfView = peel.camera
        view.allowsCameraControl = true
        view.defaultCameraController.interactionMode = .orbitTurntable
        view.defaultCameraController.worldUp = SCNVector3(0, 1, 0)
        view.defaultCameraController.target = peel.target
        view.antialiasingMode = .multisampling4X
        view.isPlaying = true   // keep rendering, so the fades animate
        window.contentView!.addSubview(view)

        let info = label(15, NSColor(white: 0.95, alpha: 1)), keys = label(12, NSColor(srgbRed: 0.65, green: 0.68, blue: 0.72, alpha: 1))
        keys.stringValue = "drag rotate   scroll/pinch zoom   |   space play/pause   left/right step"
            + "   up/down speed   0 restart   c start view   r fit   f fly to cursor   q quit"
        for t in [info, keys] { window.contentView!.addSubview(t) }
        NSLayoutConstraint.activate([
            info.leadingAnchor.constraint(equalTo: window.contentView!.leadingAnchor, constant: 18),
            info.topAnchor.constraint(equalTo: window.contentView!.topAnchor, constant: 14),
            keys.leadingAnchor.constraint(equalTo: window.contentView!.leadingAnchor, constant: 18),
            keys.bottomAnchor.constraint(equalTo: window.contentView!.bottomAnchor, constant: -12),
        ])
        peeler.onChange = { info.stringValue = self.peeler.text }
        peeler.onChange()

        let aim = { view.defaultCameraController.target = self.peel.target }
        view.keys = [
            "space": peeler.toggle, "Right": peeler.forward, "Left": peeler.back,
            "Up": { self.peeler.speed(1 / 1.5) }, "Down": { self.peeler.speed(1.5) }, "0": peeler.restart,
            "c": { self.peel.startView(); aim() },
            "r": { self.peel.frame(nodes: self.peel.layerNodes); aim() },
            "f": { self.flyToCursor(view) },
            "q": { NSApp.terminate(nil) },
        ]
        timer = Timer.scheduledTimer(withTimeInterval: 0.05, repeats: true) { _ in self.peeler.tick() }
        window.makeKeyAndOrderFront(nil)
        window.makeFirstResponder(view)
        NSApp.activate(ignoringOtherApps: true)
        print(String(format: "%.2f M triangles in %d layers, ready in %.1f s", Double(peel.triangles) / 1e6,
                     peel.layerNodes.count, Date().timeIntervalSince(t0)))
    }

    // f: centre the view on the surface under the cursor and move a third of the way to it
    func flyToCursor(_ view: SCNView) {
        let p = view.convert(window.mouseLocationOutsideOfEventStream, from: nil)
        guard let hit = view.hitTest(p, options: [.searchMode: SCNHitTestSearchMode.closest.rawValue, .ignoreHiddenNodes: true]).first,
              let cam = view.pointOfView else { return }
        let w = hit.worldCoordinates, c = cam.worldPosition
        SCNTransaction.begin()
        SCNTransaction.animationDuration = 0.6
        cam.worldPosition = SCNVector3(c.x + (w.x - c.x) / 3, c.y + (w.y - c.y) / 3, c.z + (w.z - c.z) / 3)
        cam.look(at: w, up: SCNVector3(0, 1, 0), localFront: SCNVector3(0, 0, -1))
        SCNTransaction.commit()
        view.defaultCameraController.target = w
        peel.target = w
    }

    // --snapshot out.png [--peel N] [--size WxH]: render one frame offscreen and quit
    func snapshot(to path: String, _ args: [String]) {
        var peelN = 0, size = CGSize(width: 1500, height: 1000)
        if let i = args.firstIndex(of: "--peel"), i + 1 < args.count { peelN = Int(args[i + 1]) ?? 0 }
        if let i = args.firstIndex(of: "--size"), i + 1 < args.count {
            let wh = args[i + 1].split(separator: "x").compactMap { Double($0) }
            if wh.count == 2 { size = CGSize(width: wh[0], height: wh[1]) }
        }
        for i in 0..<min(peelN, peel.layerNodes.count - 1) { peel.layerNodes[i].isHidden = true }
        if peelN > 0 { peel.frame(nodes: peel.layerNodes) }   // as r does
        let r = SCNRenderer(device: MTLCreateSystemDefaultDevice(), options: nil)
        r.scene = peel.scene
        r.pointOfView = peel.camera
        let image = r.snapshot(atTime: 0, with: size, antialiasingMode: .multisampling4X)
        if let tiff = image.tiffRepresentation, let rep = NSBitmapImageRep(data: tiff), let png = rep.representation(using: .png, properties: [:]) {
            try? png.write(to: URL(fileURLWithPath: path))
            print("wrote \(path): \(peel.triangles) triangles, \(peelN) layers peeled")
        }
        NSApp.terminate(nil)
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ app: NSApplication) -> Bool { true }
}

@main
enum AutopeelApp {
    static func main() {
        let app = NSApplication.shared
        let delegate = AppDelegate()
        app.delegate = delegate
        app.setActivationPolicy(.regular)
        let menu = NSMenu(), item = NSMenuItem()
        menu.addItem(item)
        let appMenu = NSMenu()
        let name = ProcessInfo.processInfo.processName
        appMenu.addItem(withTitle: "Quit \(name)", action: #selector(NSApplication.terminate(_:)), keyEquivalent: "q")
        item.submenu = appMenu
        app.mainMenu = menu
        app.run()
    }
}
