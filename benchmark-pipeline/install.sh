#!/usr/bin/env bash
# One-shot setup for the riker benchmark pipeline.
#
# Idempotent: re-running is safe and only re-does steps that are out of date.
# Sets up an isolated environment under benchmark-pipeline/.rust/ so the
# host's ~/.cargo and ~/.rustup are never touched.
#
#   ./install.sh                  # install everything (default)
#   ./install.sh --system-rust    # use whatever cargo is on PATH instead of
#                                 # installing rustup into .rust/
#   ./install.sh --skip-build     # set up envs but don't (re)build riker
#
# After this script succeeds, run benchmarks with `./run.sh`.
set -euo pipefail

cd "$(dirname "$0")"
PIPELINE_DIR="$PWD"
REPO_ROOT="$(cd .. && pwd)"

USE_ISOLATED_RUST=1
DO_BUILD=1
for arg in "$@"; do
  case "$arg" in
    --system-rust) USE_ISOLATED_RUST=0 ;;
    --skip-build)  DO_BUILD=0 ;;
    -h|--help)
      sed -n '2,15p' "$0"
      exit 0
      ;;
    *)
      echo "unknown flag: $arg" >&2
      exit 2
      ;;
  esac
done

log() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }

# ---- Build toolchain ------------------------------------------------------
# Slim Linux AMIs (AL2023 minimal etc.) ship without a C/C++ toolchain,
# cmake, or libclang. We need:
#   - gcc / g++ / cmake: rustc build scripts + Rust crates with C deps
#     (e.g. libz-ng-sys -> cmake)
#   - libclang: bindgen (used by hts-sys) dlopens libclang.so at build time
# Install all up front when missing so the rust + cargo-multivers + riker
# chain works on a fresh box.
have_cc=true
have_cxx=true
have_cmake=true
have_libclang=true
command -v cc    >/dev/null 2>&1 || command -v gcc >/dev/null 2>&1 || have_cc=false
command -v c++   >/dev/null 2>&1 || command -v g++ >/dev/null 2>&1 || have_cxx=false
command -v cmake >/dev/null 2>&1 || have_cmake=false
# libclang has no canonical command — probe the dynamic linker cache and
# common install paths.
if ! ldconfig -p 2>/dev/null | grep -q '/libclang\.so'; then
  if ! ls /usr/lib*/libclang*.so* /usr/lib/llvm*/lib/libclang*.so* 2>/dev/null | head -1 | grep -q .; then
    have_libclang=false
  fi
fi

if ! $have_cc || ! $have_cxx || ! $have_cmake || ! $have_libclang; then
  if [[ "$(uname -s)" == "Linux" ]]; then
    pkgs=()
    $have_cc    || pkgs+=(gcc)
    if command -v dnf >/dev/null 2>&1; then
      $have_cxx     || pkgs+=(gcc-c++)
      $have_cmake   || pkgs+=(cmake)
      $have_libclang || pkgs+=(clang-devel)
      log "Installing build toolchain: ${pkgs[*]}"
      sudo dnf install -y -q "${pkgs[@]}"
    elif command -v apt-get >/dev/null 2>&1; then
      # When c++ is missing we install build-essential, which covers BOTH
      # gcc and g++ — drop any standalone gcc from the list to avoid
      # apt-get complaining about a redundant package.
      $have_cxx || pkgs=("${pkgs[@]/gcc/}")
      $have_cxx || pkgs+=(build-essential)
      $have_cmake   || pkgs+=(cmake)
      $have_libclang || pkgs+=(libclang-dev)
      log "Installing build toolchain: ${pkgs[*]}"
      sudo apt-get update -qq
      sudo apt-get install -y -qq "${pkgs[@]}"
    else
      echo "ERROR: missing build toolchain and unrecognized package manager." >&2
      echo "Install gcc, g++, cmake, and libclang manually and re-run." >&2
      exit 1
    fi
  else
    echo "ERROR: missing build toolchain on a non-Linux host." >&2
    echo "On macOS: \`xcode-select --install && brew install cmake llvm\`." >&2
    exit 1
  fi
fi

# ---- pixi -----------------------------------------------------------------
if ! command -v pixi >/dev/null 2>&1; then
  log "Installing pixi to ~/.pixi"
  curl -fsSL https://pixi.sh/install.sh | bash
  export PATH="$HOME/.pixi/bin:$PATH"
else
  log "pixi found: $(pixi --version)"
fi

# ---- pixi environment -----------------------------------------------------
log "Materializing pixi environment (single 'default' env)"
pixi install

# ---- rust toolchain -------------------------------------------------------
RUST_DIR="$PIPELINE_DIR/.rust"
if [[ "$USE_ISOLATED_RUST" -eq 1 ]]; then
  export RUSTUP_HOME="$RUST_DIR/rustup"
  export CARGO_HOME="$RUST_DIR/cargo"
  export PATH="$CARGO_HOME/bin:$PATH"

  if [[ ! -x "$CARGO_HOME/bin/rustup" ]]; then
    log "Installing rustup into $RUST_DIR (host ~/.cargo and ~/.rustup are not touched)"
    mkdir -p "$RUST_DIR"
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
      | sh -s -- -y --no-modify-path --default-toolchain none --profile minimal
  else
    log "Isolated rustup present at $CARGO_HOME"
  fi

  # rust-toolchain.toml at the repo root pins the channel; force its install
  # now so build progress is visible up front rather than mixed into the
  # riker build output.
  ( cd "$REPO_ROOT" && rustup show active-toolchain >/dev/null )
else
  if ! command -v cargo >/dev/null 2>&1; then
    echo "ERROR: --system-rust set but cargo not on PATH. Install rustup or drop --system-rust." >&2
    exit 1
  fi
  log "Using system cargo: $(cargo --version)"
fi

# ---- cargo-multivers (x86_64 only) ----------------------------------------
ARCH="$(uname -m)"
if [[ "$ARCH" == "x86_64" || "$ARCH" == "amd64" ]]; then
  if ! command -v cargo-multivers >/dev/null 2>&1; then
    log "Installing cargo-multivers"
    # Floor of 0.12.0: carries the fexecve/memfd fix needed for the launcher
    # to run under binfmt emulation (fulcrumgenomics/riker#32).
    cargo install cargo-multivers --locked --version '>=0.12.0'
  else
    log "cargo-multivers found"
  fi
fi

# ---- build riker ----------------------------------------------------------
# x86_64: cargo-multivers builds v1/v2/v4 variants and a launcher. This
# matches the artifact users actually `cargo install`, so it's the right
# binary to benchmark.
# aarch64: cargo build --release. cargo-multivers' default CPU list is
# x86_64-only and the Cargo.toml [package.metadata.multivers.x86_64] block
# limits it to that; on aarch64 generic codegen is what users get.
if [[ "$DO_BUILD" -eq 1 ]]; then
  cd "$REPO_ROOT"
  case "$ARCH" in
    x86_64|amd64)
      log "Building riker with cargo-multivers (x86-64 v1/v2/v4 variants)"
      cargo multivers --profile dist
      # cargo-multivers writes the launcher under target/<target-triple>/<profile>/.
      # Symlink to a stable path so config.yaml can use ../target/dist/riker
      # regardless of host triple.
      MV_LAUNCHER="$(find target -type f -name riker \
                       -not -path '*/deps/*' -not -path '*/build/*' \
                       -path '*/dist/riker' \
                       -not -path 'target/dist/riker' \
                       -print -quit)"
      if [[ -z "${MV_LAUNCHER:-}" ]]; then
        echo "ERROR: cargo multivers reported success but no launcher binary was found under target/" >&2
        exit 1
      fi
      mkdir -p target/dist
      ln -sf "../../$MV_LAUNCHER" target/dist/riker
      ;;
    aarch64|arm64)
      log "Building riker with cargo build --release (single aarch64 binary)"
      cargo build --release
      mkdir -p target/dist
      ln -sf "../release/riker" target/dist/riker
      ;;
    *)
      echo "ERROR: unsupported architecture: $ARCH" >&2
      exit 1
      ;;
  esac
  cd "$PIPELINE_DIR"
  log "riker binary: $(realpath ../target/dist/riker)"
  ../target/dist/riker --version || true
fi

cat <<EOF

Setup complete.

Next:
  1. Smoke-test (~30 min, downloads the 4x downsample only):
       ./run.sh config/smoke.config.yaml

  2. Full performance run (heavy: ~180 GB of staged BAMs, several hours):
       ./run.sh config/performance.config.yaml

  3. Preview the job graph without running:
       ./run.sh --dry-run

EOF
