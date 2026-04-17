#!/usr/bin/env bash
# Build the upstream lofreq C binary from source for parity testing.
#
# Pin:    v2.1.5  (last tagged release at time of writing)
# Output: parity/upstream/bin/lofreq
#
# Dependencies (apt-installable on Debian/Ubuntu):
#   build-essential autoconf zlib1g-dev libbz2-dev liblzma-dev
#   libcurl4-openssl-dev libssl-dev libncurses5-dev wget git
#
# The script is idempotent: if parity/upstream/bin/lofreq already exists
# and --force is not set, it exits early.

set -euo pipefail

LOFREQ_VERSION="${LOFREQ_VERSION:-2.1.5}"
HTSLIB_VERSION="${HTSLIB_VERSION:-1.20}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PARITY="$ROOT/parity/upstream"
BIN="$PARITY/bin/lofreq"
BUILD="$PARITY/build"

force=0
if [[ "${1:-}" == "--force" ]]; then force=1; fi

if [[ -x "$BIN" && $force -eq 0 ]]; then
  echo "[build-upstream] $BIN already present. Pass --force to rebuild."
  "$BIN" version | head -2 || true
  exit 0
fi

mkdir -p "$PARITY/bin" "$BUILD"
cd "$BUILD"

# --- htslib ------------------------------------------------------------
if [[ ! -d "htslib-$HTSLIB_VERSION" ]]; then
  echo "[build-upstream] fetching htslib $HTSLIB_VERSION"
  wget -q "https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2"
  tar -xjf "htslib-$HTSLIB_VERSION.tar.bz2"
fi
pushd "htslib-$HTSLIB_VERSION" >/dev/null
if [[ ! -f libhts.a ]]; then
  ./configure --prefix="$PARITY" --disable-libcurl --disable-bz2 --disable-lzma
  make -j"$(nproc)"
  make install
fi
popd >/dev/null

# --- lofreq ------------------------------------------------------------
if [[ ! -d "lofreq-$LOFREQ_VERSION" ]]; then
  echo "[build-upstream] fetching lofreq $LOFREQ_VERSION"
  wget -q "https://github.com/CSB5/lofreq/archive/refs/tags/v$LOFREQ_VERSION.tar.gz" \
    -O "lofreq-$LOFREQ_VERSION.tar.gz"
  tar -xzf "lofreq-$LOFREQ_VERSION.tar.gz"
fi
cd "lofreq-$LOFREQ_VERSION"
if [[ ! -x src/lofreq/lofreq ]]; then
  ./bootstrap
  ./configure \
    --prefix="$PARITY" \
    --with-htslib="$PARITY" \
    SAMTOOLS="" BCFTOOLS=""
  make -j"$(nproc)"
fi

# `make install` in the lofreq tree also tries to install its Python
# helper scripts via setuptools, which fails on newer setuptools versions.
# We only need the C binary; copy it into place manually.
cp "src/lofreq/lofreq" "$BIN"

echo "[build-upstream] installed:"
"$BIN" version | head -2
