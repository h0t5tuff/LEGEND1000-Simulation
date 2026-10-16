// LEGEND.app: one app for the repo's three viewers. On launch it offers l200 and l1000 (autopeel.swift: the LEGEND-200
// and LEGEND-1000 geometries, peeled layer by layer) and rt (rt.swift: the LEGEND-1000 re-entrance tube designs);
// each reads its own folder of the bundle's Resources, all written by make_app.py. Closing a viewer comes back here.
// Built by make_app.py:
//   xcrun swiftc -O -parse-as-library -swift-version 5 legend.swift autopeel.swift rt.swift
// Hidden check: <app>/Contents/MacOS/LEGEND --snapshot out.png [--open l200|l1000|rt] [that viewer's options] draws
// the choices (or that viewer) offscreen and quits.

import AppKit
import SceneKit

protocol Viewer: AnyObject {
    var window: NSWindow { get }
    func show()
    func close()   // stop its timers and keys, so the viewer and its scene can go
    func snapshot(to path: String, _ args: [String])
}

func argument(_ args: [String], _ flag: String) -> String? {
    args.firstIndex(of: flag).flatMap { $0 + 1 < args.count ? args[$0 + 1] : nil }
}

func parseSize(_ s: String) -> NSSize? {
    let wh = s.split(separator: "x").compactMap { Double($0) }
    return wh.count == 2 ? NSSize(width: wh[0], height: wh[1]) : nil
}

func writePNG(_ image: NSImage, _ path: String) {
    if let tiff = image.tiffRepresentation, let rep = NSBitmapImageRep(data: tiff),
       let png = rep.representation(using: .png, properties: [:]) {
        try? png.write(to: URL(fileURLWithPath: path))
    }
}

// ── shared by the viewers and the choices ─────────────────────────────────────────────────────────────

let fg = NSColor(srgbRed: 0.91, green: 0.92, blue: 0.93, alpha: 1)
let muted = NSColor(srgbRed: 0.60, green: 0.64, blue: 0.70, alpha: 1)
let panelColour = NSColor(srgbRed: 0.08, green: 0.09, blue: 0.11, alpha: 1)
let cardColour = NSColor(srgbRed: 0.11, green: 0.125, blue: 0.153, alpha: 1)
let lineColour = NSColor(srgbRed: 0.20, green: 0.227, blue: 0.27, alpha: 1)

func gradient(bottom: NSColor, top: NSColor) -> NSImage {
    let image = NSImage(size: NSSize(width: 8, height: 512))
    image.lockFocus()
    NSGradient(starting: bottom, ending: top)!.draw(in: NSRect(x: 0, y: 0, width: 8, height: 512), angle: 90)
    image.unlockFocus()
    return image
}

final class KeyView: SCNView {   // a 3D view whose keys call functions
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

class Fill: NSView {   // a panel's background
    override func draw(_ dirty: NSRect) { panelColour.setFill(); bounds.fill() }
}

class Card: NSView {   // a rounded box
    override func draw(_ dirty: NSRect) {
        let p = NSBezierPath(roundedRect: bounds.insetBy(dx: 0.5, dy: 0.5), xRadius: 10, yRadius: 10)
        cardColour.setFill(); p.fill()
        lineColour.setStroke(); p.lineWidth = 1; p.stroke()
    }
}

func field(_ size: CGFloat, _ weight: NSFont.Weight = .regular, _ colour: NSColor = fg, _ s: String = "") -> NSTextField {
    let t = NSTextField(labelWithString: s)
    t.font = .monospacedDigitSystemFont(ofSize: size, weight: weight)
    t.textColor = colour
    return t
}

func wrapping(_ s: String, _ size: CGFloat, _ colour: NSColor, width: CGFloat) -> NSTextField {
    let t = NSTextField(wrappingLabelWithString: s)
    t.font = .systemFont(ofSize: size)
    t.textColor = colour
    t.preferredMaxLayoutWidth = width
    t.isSelectable = false
    return t
}

func row(_ a: NSView, _ b: NSView) -> NSStackView {   // a at the left edge, b at the right
    let s = NSStackView(views: [a, b])
    s.orientation = .horizontal
    s.distribution = .equalSpacing
    return s
}

func column(_ views: [NSView], spacing: CGFloat, width: CGFloat) -> NSStackView {
    let s = NSStackView(views: views)
    s.orientation = .vertical
    s.alignment = .leading
    s.spacing = spacing
    s.translatesAutoresizingMaskIntoConstraints = false
    for v in views {
        v.translatesAutoresizingMaskIntoConstraints = false
        v.widthAnchor.constraint(equalToConstant: width).isActive = true
    }
    return s
}

func card(_ views: [NSView], width: CGFloat) -> NSView {
    let c = Card(), s = column(views, spacing: 8, width: width - 32)
    c.addSubview(s)
    NSLayoutConstraint.activate([
        s.leadingAnchor.constraint(equalTo: c.leadingAnchor, constant: 16),
        s.topAnchor.constraint(equalTo: c.topAnchor, constant: 14),
        s.bottomAnchor.constraint(equalTo: c.bottomAnchor, constant: -14),
    ])
    return c
}

// ── the choices ─────────────────────────────────────────────────────────────────────────────────────

let choices: [(key: String, line: String)] = [
    ("l200", "The LEGEND-200 geometry, peeled layer by layer as you rotate it."),
    ("l1000", "The LEGEND-1000 geometry, peeled layer by layer as you rotate it."),
    ("rt", "The LEGEND-1000 re-entrance tube: the BI of every split into SS, Cu and EFCu."),
]

final class ChoiceCard: Card {
    var pick: () -> Void = {}
    override func mouseUp(with event: NSEvent) { if bounds.contains(convert(event.locationInWindow, from: nil)) { pick() } }
    override func resetCursorRects() { addCursorRect(bounds, cursor: .pointingHand) }
}

final class ChooserView: Fill {
    var keys: [String: () -> Void] = [:]
    override var acceptsFirstResponder: Bool { true }
    override func keyDown(with event: NSEvent) {
        if let f = keys[event.charactersIgnoringModifiers ?? ""] { f() } else { super.keyDown(with: event) }
    }
}

final class AppDelegate: NSObject, NSApplicationDelegate, NSWindowDelegate {
    var chooser: NSWindow!
    let status = field(13, .regular, muted)
    var viewer: Viewer?

    func applicationDidFinishLaunching(_ note: Notification) {
        makeChooser()
        let args = CommandLine.arguments
        if let path = argument(args, "--snapshot") {
            if let key = argument(args, "--open") {
                if choices.contains(where: { $0.key == key }), let v = make(key) { v.snapshot(to: path, args) }
                else { print("--open takes one of: \(choices.map { $0.key }.joined(separator: ", "))") }
            } else {
                let root = chooser.contentView!
                root.layoutSubtreeIfNeeded()
                if let rep = root.bitmapImageRepForCachingDisplay(in: root.bounds) {
                    root.cacheDisplay(in: root.bounds, to: rep)
                    let image = NSImage(size: root.bounds.size)
                    image.addRepresentation(rep)
                    writePNG(image, path)
                    print("wrote \(path): the choices")
                }
            }
            NSApp.terminate(nil)
            return
        }
        chooser.center()
        chooser.makeKeyAndOrderFront(nil)
        chooser.makeFirstResponder(chooser.contentView)
        NSApp.activate(ignoringOtherApps: true)
    }

    func makeChooser() {
        let cardW: CGFloat = 250, cardH: CGFloat = 372, gap: CGFloat = 20, pad: CGFloat = 32
        let width = 3 * cardW + 2 * gap, size = NSSize(width: width + 2 * pad, height: cardH + 150)
        chooser = NSWindow(contentRect: NSRect(origin: .zero, size: size), styleMask: [.titled, .closable, .miniaturizable],
                           backing: .buffered, defer: false)
        chooser.isReleasedWhenClosed = false
        chooser.title = "LEGEND"
        chooser.appearance = NSAppearance(named: .darkAqua)
        chooser.delegate = self
        let root = ChooserView(frame: NSRect(origin: .zero, size: size))
        chooser.contentView = root
        let res = Bundle.main.resourceURL!
        let cards: [NSView] = choices.enumerated().map { k, c in
            let card = ChoiceCard()
            let image = NSImageView(image: NSImage(contentsOf: res.appendingPathComponent("\(c.key)/render.png")) ?? NSImage())
            image.imageScaling = .scaleProportionallyUpOrDown
            image.heightAnchor.constraint(equalToConstant: 210).isActive = true
            let s = column([image, field(22, .bold, fg, c.key), wrapping(c.line, 13, muted, width: cardW - 32),
                            field(11, .regular, muted, "press \(k + 1)")], spacing: 8, width: cardW - 32)
            card.addSubview(s)
            card.translatesAutoresizingMaskIntoConstraints = false
            NSLayoutConstraint.activate([
                card.widthAnchor.constraint(equalToConstant: cardW), card.heightAnchor.constraint(equalToConstant: cardH),
                s.leadingAnchor.constraint(equalTo: card.leadingAnchor, constant: 16),
                s.topAnchor.constraint(equalTo: card.topAnchor, constant: 16),
            ])
            card.pick = { [unowned self] in self.open(c.key) }
            return card
        }
        let row = NSStackView(views: cards)
        row.orientation = .horizontal
        row.spacing = gap
        let all = column([field(26, .bold, fg, "LEGEND"),
                          wrapping("Click a viewer or press its number. Closing a viewer brings you back here.", 13, muted, width: width),
                          row, status], spacing: 12, width: width)
        root.addSubview(all)
        NSLayoutConstraint.activate([
            all.leadingAnchor.constraint(equalTo: root.leadingAnchor, constant: pad),
            all.topAnchor.constraint(equalTo: root.topAnchor, constant: 24),
        ])
        root.keys = ["q": { NSApp.terminate(nil) }]
        for (k, c) in choices.enumerated() { root.keys["\(k + 1)"] = { [unowned self] in self.open(c.key) } }
    }

    func make(_ key: String) -> Viewer? {
        let dir = Bundle.main.resourceURL!.appendingPathComponent(key, isDirectory: true)
        do { return key == "rt" ? try RTController(dir) : try PeelController(dir) } catch {
            let a = NSAlert()
            a.messageText = "No \(key) data in the app"
            a.informativeText = "\(error)"
            a.runModal()
            return nil
        }
    }

    // build the viewer after the chooser has drawn "opening": l1000 takes a few seconds
    func open(_ key: String) {
        guard viewer == nil else { return }
        status.stringValue = "opening \(key) …"
        DispatchQueue.main.async {
            defer { self.status.stringValue = "" }
            guard self.viewer == nil, let v = self.make(key) else { return }
            self.viewer = v
            v.window.delegate = self
            self.chooser.orderOut(nil)
            v.show()
        }
    }

    func windowWillClose(_ note: Notification) {
        guard let w = note.object as? NSWindow else { return }
        if let v = viewer, w === v.window {
            v.close()
            viewer = nil
            DispatchQueue.main.async {
                self.chooser.makeKeyAndOrderFront(nil)
                self.chooser.makeFirstResponder(self.chooser.contentView)
            }
        } else if w === chooser {
            NSApp.terminate(nil)
        }
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ app: NSApplication) -> Bool { false }
}

@main
enum LegendApp {
    static func main() {
        let app = NSApplication.shared
        let delegate = AppDelegate()
        app.delegate = delegate
        app.setActivationPolicy(.regular)
        let menu = NSMenu()
        let appItem = NSMenuItem(), appMenu = NSMenu()
        appMenu.addItem(withTitle: "Quit LEGEND", action: #selector(NSApplication.terminate(_:)), keyEquivalent: "q")
        appItem.submenu = appMenu
        let windowItem = NSMenuItem(), windowMenu = NSMenu(title: "Window")
        windowMenu.addItem(withTitle: "Close", action: #selector(NSWindow.performClose(_:)), keyEquivalent: "w")
        windowMenu.addItem(withTitle: "Minimize", action: #selector(NSWindow.performMiniaturize(_:)), keyEquivalent: "m")
        windowItem.submenu = windowMenu
        menu.addItem(appItem)
        menu.addItem(windowItem)
        app.mainMenu = menu
        app.run()
    }
}
