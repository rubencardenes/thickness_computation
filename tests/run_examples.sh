#!/usr/bin/env bash
# Runs every example command listed in README.md and checks:
#   1. the program exits successfully
#   2. its output file is produced (and is a valid PNG, for .png outputs)
#   3. for thickness tools (-m), the reported mean thickness is within 1%
#      of a known-good reference value
#
# Reference mean-thickness values were captured by running each example
# against the current algorithm and are meant to catch regressions, not to
# encode "ground truth" thickness.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUTDIR="$ROOT/tests/output"
mkdir -p "$OUTDIR"

PASS=0
FAIL=0
FAILED_NAMES=()

is_valid_png() {
  [ -s "$1" ] || return 1
  head -c 8 "$1" | od -An -tx1 | tr -d ' \n' | grep -qi '^89504e470d0a1a0a$'
}

extract_mean() {
  grep -Eo 'Mean thickness[^=]*= *[0-9.]+' "$1" | head -1 | grep -Eo '[0-9.]+$'
}

within_pct() {
  # $1=actual $2=expected $3=percent tolerance
  awk -v a="$1" -v e="$2" -v p="$3" 'BEGIN {
    diff = a - e; if (diff < 0) diff = -diff;
    tol = e * p / 100.0; if (tol < 0) tol = -tol;
    exit !(diff <= tol)
  }'
}

fail() {
  echo "FAIL: $1"
  FAIL=$((FAIL + 1))
  FAILED_NAMES+=("$1")
}

pass() {
  echo "PASS: $1"
  PASS=$((PASS + 1))
}

# run_test <name> <command...> -- runs the command, logs it, checks exit code.
# On success prints nothing further; caller does its own output/mean checks.
run_cmd() {
  local name="$1"; shift
  local log="$OUTDIR/${name}.log"
  echo "=== $name ==="
  echo "\$ $*" > "$log"
  if ! "$@" >> "$log" 2>&1; then
    fail "$name (command exited non-zero, see $log)"
    return 1
  fi
  return 0
}

check_output_file() {
  local name="$1" outfile="$2"
  if [[ "$outfile" == *.png ]]; then
    if ! is_valid_png "$outfile"; then
      fail "$name (missing/invalid PNG output: $outfile)"
      return 1
    fi
  else
    if [ ! -s "$outfile" ]; then
      fail "$name (missing/empty output: $outfile)"
      return 1
    fi
  fi
  return 0
}

check_mean() {
  local name="$1" log="$2" expected="$3"
  local actual
  actual=$(extract_mean "$log")
  if [ -z "$actual" ]; then
    fail "$name (no 'Mean thickness' line found in $log)"
    return 1
  fi
  if ! within_pct "$actual" "$expected" 1; then
    fail "$name (mean thickness $actual not within 1% of expected $expected)"
    return 1
  fi
  return 0
}

# thickness_test <name> <outfile> <expected_mean> <command...>
thickness_test() {
  local name="$1" outfile="$2" expected="$3"; shift 3
  run_cmd "$name" "$@" || return
  check_output_file "$name" "$outfile" || return
  check_mean "$name" "$OUTDIR/${name}.log" "$expected" || return
  pass "$name"
}

# plain_test <name> <outfile> <command...> -- no mean-thickness check
plain_test() {
  local name="$1" outfile="$2"; shift 2
  run_cmd "$name" "$@" || return
  check_output_file "$name" "$outfile" || return
  pass "$name"
}

echo "Building project (make)..."
if ! make -C "$ROOT" > "$OUTDIR/build.log" 2>&1; then
  echo "BUILD FAILED, see $OUTDIR/build.log"
  cat "$OUTDIR/build.log"
  exit 1
fi
echo "Build OK."
echo

### thickness2D examples ######################################################

thickness_test thickness2d_ring1 "$OUTDIR/thickness2d_ring1.png" 7.677905 \
  ./thickness2D -n 20 -i 100 -c 2 -s -m --hx 0.5 --hy 0.5 --lw 40 --lc 98 \
  data/domain_anillo_256_1.png "$OUTDIR/thickness2d_ring1.png"

thickness_test thickness2d_ring2 "$OUTDIR/thickness2d_ring2.png" 11.702454 \
  ./thickness2D -n 20 -i 100 -c 2 -s -m --hx 0.5 --hy 0.5 --lw 40 --lc 98 \
  data/domain_anillo_256_2.png "$OUTDIR/thickness2d_ring2.png"

thickness_test thickness2d_ellipse "$OUTDIR/thickness2d_ellipse.png" 17.047857 \
  ./thickness2D -n 20 -i 100 -c 2 -s -m --hx 0.5 --hy 0.5 --lw 40 --lc 98 \
  data/domain_elipse.png "$OUTDIR/thickness2d_ellipse.png"

thickness_test thickness2d_square "$OUTDIR/thickness2d_square.png" 32.201904 \
  ./thickness2D -n 30 -i 200 -c 2 -s -m --hx 0.5 --hy 0.5 --lw 40 --lc 98 \
  data/domain_cuadrado.png "$OUTDIR/thickness2d_square.png"

thickness_test thickness2d_donut "$OUTDIR/thickness2d_donut.png" 14.532133 \
  ./thickness2D -n 15 -i 100 -c 2 -s -m --hx 0.5 --hy 0.5 --lw 40 --lc 98 \
  data/domain_donut.png "$OUTDIR/thickness2d_donut.png"

### thickness2D_knee example ##################################################

thickness_test thickness2d_knee "$OUTDIR/thickness2d_knee.png" 5.605474 \
  ./thickness2D_knee -c 2 -i 200 -n 10 --hx 1 --hy 1 -m \
  data/domain_knee_512.png "$OUTDIR/thickness2d_knee.png"

### thickness3D examples ######################################################

# thickness3d_caja and thickness3d_elipsoid reference values below were
# recaptured after fixing a real i/j-axis transposition bug in
# laplace3D_voxelsize/laplace3D (they read the domain buffer with a
# different stride convention than compute_boundary_cortex3D/EdgeDetect3D
# use to populate it). Both phantoms have i/j-asymmetric geometry (the caja
# core spans i in [30,49] but j in [35,44]; the ellipsoid has rx=15 vs
# ry=35), so the fix changes their measured thickness; thickness3d_sphere is
# isotropic and is unaffected, confirming the change is a correctness fix
# and not a regression.
#
# thickness3d_elipsoid's value was recaptured again after pinning
# -ffp-contract=off in the Makefile: laplace3D_voxelsize's relaxation runs a
# fixed 200 iterations without converging, so letting the compiler fuse
# a*b+c into one rounding step made the result depend on whether the build
# target has hardware FMA (arm64 always does; x86_64 without -mfma usually
# doesn't), moving this anisotropic phantom's mean by >1% between a macOS
# and a GitHub Actions (ubuntu-latest/x86_64) build of identical source.
# thickness3d_caja's value happened not to move.
thickness_test thickness3d_caja "$OUTDIR/thickness3d_caja.volf" 23.825426 \
  ./thickness3D -m -n 20 -i 200 --lw 3 --lc 2 \
  data/input_caja3d.vols "$OUTDIR/thickness3d_caja.volf" 80 80 80

thickness_test thickness3d_elipsoid "$OUTDIR/thickness3d_elipsoid.volf" 12.575148 \
  ./thickness3D -m -n 20 -i 200 --lw 3 --lc 2 \
  data/phantom_elipsoid.vols "$OUTDIR/thickness3d_elipsoid.volf" 80 80 80

thickness_test thickness3d_sphere "$OUTDIR/thickness3d_sphere.volf" 24.552910 \
  ./thickness3D -m -n 20 -i 200 --lw 3 --lc 2 \
  data/phantom_sphere.vols "$OUTDIR/thickness3d_sphere.volf" 80 80 80

### laplace2D examples (no mean-thickness output) #############################

plain_test laplace2d_png "$OUTDIR/laplace2d.png" \
  ./laplace2D data/domain_anillo_poisson.png "$OUTDIR/laplace2d.png" 100

plain_test laplace2d_flt "$OUTDIR/laplace2d.flt" \
  ./laplace2D data/domain_anillo_poisson.png "$OUTDIR/laplace2d.flt" 100

### laplace3D example ##########################################################

plain_test laplace3d "$OUTDIR/laplace3d.volf" \
  ./laplace3D 80 80 80 data/domain_anillo_3d.vol "$OUTDIR/laplace3d.volf" 100

### poisson2D example ##########################################################

plain_test poisson2d "$OUTDIR/poisson2d.png" \
  ./poisson2D data/domain_anillo_poisson.png "$OUTDIR/poisson2d.png" 200

# poisson2D always additionally writes max_local.png in the cwd (the project
# root, since we run from there) -- check that side-effect output too.
if check_output_file poisson2d_max_local "$ROOT/max_local.png"; then
  pass poisson2d_max_local
fi

### unit checks ###############################################################

# Connectivity rules for the flood fills (see tests/check_adjacency.c).
echo "=== check_adjacency ==="
if make -C "$ROOT" check_adjacency > "$OUTDIR/check_adjacency.log" 2>&1 \
   && "$ROOT/check_adjacency" >> "$OUTDIR/check_adjacency.log" 2>&1; then
  pass check_adjacency
else
  cat "$OUTDIR/check_adjacency.log"
  fail check_adjacency
fi

################################################################################

echo
echo "==================================="
echo "Passed: $PASS  Failed: $FAIL"
if [ $FAIL -gt 0 ]; then
  echo "Failed tests: ${FAILED_NAMES[*]}"
  exit 1
fi
exit 0
