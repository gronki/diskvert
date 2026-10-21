#!/usr/bin/env bash
# Builds a diskvert AppImage (the diskvert program only). Runs inside Ubuntu 20.04,
# so the AppImage works on systems with glibc 2.31 or newer.
#
# Local usage (Docker), from the repository root:
#   docker run --rm -v $PWD:/src ubuntu:20.04 bash /src/scripts/build_appimage.sh
#
# Output: diskvert-<arch>.AppImage in the repository root (or $OUT_DIR).

set -euo pipefail

SRC_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="${OUT_DIR:-${SRC_DIR}}"

export DEBIAN_FRONTEND=noninteractive

apt-get update -y
apt-get install -y --no-install-recommends \
    make gfortran libopenblas-dev wget file ca-certificates

# ── build (in a copy, so the source tree is not touched) ──────────────────
BUILD_DIR="$(mktemp -d)"
cp -a "$SRC_DIR/." "$BUILD_DIR/"
make -C "$BUILD_DIR/build" clean
make -C "$BUILD_DIR/build" diskvert

APPDIR="$(mktemp -d)/AppDir"
mkdir -p "$APPDIR/usr/bin" "$APPDIR/usr/lib"
cp "$BUILD_DIR/build/diskvert" "$APPDIR/usr/bin/"

# ── bundle shared libraries ────────────────────────────────────────────────
ldd "$APPDIR/usr/bin/diskvert" \
    | awk '/=> \// { print $3 }' \
    | grep -E 'libgfortran|libquadmath|libopenblas' \
    | sort -u \
    | xargs -I{} cp -vL {} "$APPDIR/usr/lib/"

# ── AppRun, .desktop and icon (required by appimagetool) ──────────────────
cat > "$APPDIR/AppRun" <<'EOF'
#!/usr/bin/env bash
HERE="$(cd "$(dirname "$0")" && pwd)"
export LD_LIBRARY_PATH="$HERE/usr/lib:${LD_LIBRARY_PATH:-}"
exec "$HERE/usr/bin/diskvert" "$@"
EOF
chmod +x "$APPDIR/AppRun"

cat > "$APPDIR/diskvert.desktop" <<'EOF'
[Desktop Entry]
Name=diskvert
Exec=diskvert
Icon=diskvert
Type=Application
Categories=Science;
Terminal=true
EOF

cat > "$APPDIR/diskvert.svg" <<'EOF'
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64"><rect width="64" height="64" fill="#222"/><ellipse cx="32" cy="32" rx="28" ry="8" fill="#e94626"/><circle cx="32" cy="32" r="5" fill="#000"/></svg>
EOF

# ── package ────────────────────────────────────────────────────────────────
ARCH=$(uname -m)
export ARCH
TOOL_DIR="$(mktemp -d)"
wget -q -O "$TOOL_DIR/appimagetool" \
    "https://github.com/AppImage/appimagetool/releases/download/continuous/appimagetool-${ARCH}.AppImage"
chmod +x "$TOOL_DIR/appimagetool"
cd "$TOOL_DIR" && ./appimagetool --appimage-extract > /dev/null

"$TOOL_DIR/squashfs-root/AppRun" "$APPDIR" "$OUT_DIR/diskvert-${ARCH}.AppImage"

echo ""
echo "Done: $OUT_DIR/diskvert-${ARCH}.AppImage"
