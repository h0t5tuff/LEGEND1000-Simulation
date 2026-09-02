#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# view the reentrant-tube's three sections (EFCu / OFHC / SS) live.
#
# The OFHC & SS shells share their OUTER surface with the EFCu mother (all r=931),
# which z-fights in the viewer. This makes a THROWAWAY copy of l1000.gdml with those
# two shells nudged +2 mm outward (so they win the depth test), opens the interactive
# viewer isolated + framed on the RT, then deletes the copy on exit. l1000.gdml is
# never touched, and nothing 18 MB is left in the repo.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

command -v remage >/dev/null 2>&1 || {
  echo "error: 'remage' not on PATH. Open your normal terminal (where 'remage --version'" >&2
  echo "       works) and run ./l1000-rt.sh from there." >&2
  exit 1
}
[ -r "$DIR/l1000.gdml" ]       || { echo "error: l1000.gdml not found beside this script." >&2; exit 1; }
[ -r "$DIR/l1000-colors.mac" ] || { echo "error: l1000-colors.mac not found beside this script." >&2; exit 1; }

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
GDML="$TMP/l1000-rt.gdml"
MAC="$TMP/l1000-rt-view.mac"

# nudge ONLY the two daughter outer-bound radii (930.999989998 -> 933.0; every one of
# the 85 matches lives inside ofhc_cu_outer_bound + ss_316l_outer_bound — surgical).
sed 's/930\.999989998/933.0/g' "$DIR/l1000.gdml" > "$GDML"

# geometry via IncludeGDMLFile — NOT `remage -g`: that flag is variadic
# (--gdml-files FILE ...) and would swallow the macro arg, so nothing initializes.
cat > "$MAC" <<EOF
/RMG/Geometry/GDMLDisableOverlapCheck true
/RMG/Geometry/IncludeGDMLFile $GDML
/run/initialize
EOF
cat >> "$MAC" <<'EOF'
/vis/open OGL
/vis/viewer/set/lineSegmentsPerCircle 24
/control/execute l1000-colors.mac
/vis/drawVolume
/vis/viewer/set/style surface
/vis/viewer/set/background 0 0 0
# isolate the RT: hide everything, show only the three sections
/vis/geometry/set/visibility world -1 false
/vis/geometry/set/visibility reentrancetube 0 true
/vis/geometry/set/visibility ofhc_cu 0 true
/vis/geometry/set/visibility ss_316l 0 true
# frame it upright, side-on, centred on the RT (z-centre ~1.9 m)
/vis/viewer/set/upVector 0 0 1
/vis/viewer/set/viewpointThetaPhi 90 0 deg
/vis/viewer/set/lightsThetaPhi 70 15 deg
/vis/viewer/set/targetPoint 0 0 1.9 m
/vis/viewer/zoom 0.8
/vis/viewer/set/autoRefresh true
EOF

echo "RT viewer — EFCu (copper) / OFHC (bronze) / SS (blue-grey); seams at z = 852 & 2857 mm."
echo "  rotate/zoom freely.   bring the full geometry back:  /vis/geometry/set/visibility world -1 true"
remage -i "$MAC"
