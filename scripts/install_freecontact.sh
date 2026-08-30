#!/usr/bin/env bash
# Install freecontact into the active conda environment, without root.
#
# freecontact is not packaged for conda on any channel (rate4site is: `conda install -c hcc
# rate4site`, and it is already in environment.yml). The only upstream binary distribution is the
# Debian/Ubuntu package, so this fetches that and unpacks it into $CONDA_PREFIX.
#
# If you have root, `sudo apt-get install freecontact` is simpler and does the same thing.
# This script exists for machines and CI runners where that is not available.
#
# Usage:
#   conda activate thoipapy
#   ./scripts/install_freecontact.sh
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: no conda environment is active. Run 'conda activate thoipapy' first." >&2
    exit 1
fi

PACKAGES=(freecontact libfreecontact0t64 libboost-program-options1.83.0 libxerces-c3.2t64)
workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

echo "Downloading: ${PACKAGES[*]}"
(cd "$workdir" && apt-get download "${PACKAGES[@]}")

echo "Unpacking into $CONDA_PREFIX"
for deb in "$workdir"/*.deb; do
    dpkg-deb -x "$deb" "$workdir/root"
done
mkdir -p "$CONDA_PREFIX/libexec"
cp -a "$workdir/root/usr/bin/freecontact" "$CONDA_PREFIX/libexec/freecontact.bin"
cp -a "$workdir/root/usr/lib/x86_64-linux-gnu/"*.so* "$CONDA_PREFIX/lib/"

# The Debian binary has no RPATH, and the dynamic loader does not search $CONDA_PREFIX/lib for
# binaries that conda did not build. Rather than depend on patchelf, install a wrapper on PATH
# that points the loader at the environment's own lib directory.
cat > "$CONDA_PREFIX/bin/freecontact" <<'WRAPPER'
#!/usr/bin/env bash
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export LD_LIBRARY_PATH="$here/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$here/libexec/freecontact.bin" "$@"
WRAPPER
chmod +x "$CONDA_PREFIX/bin/freecontact"

if ! "$CONDA_PREFIX/bin/freecontact" --version >/dev/null 2>&1; then
    echo "ERROR: freecontact was installed but does not run." >&2
    exit 1
fi
echo "freecontact $("$CONDA_PREFIX/bin/freecontact" --version) installed to $CONDA_PREFIX/bin"
