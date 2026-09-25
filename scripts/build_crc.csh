#!/bin/tcsh
#
# Build AZURE2 (GUI + CLI + pyazr) on the Notre Dame CRC cluster.
#
# Usage:
#   cd ~/AZURE2          (wherever you cloned https://github.com/AZURE2DEV/AZURE2.git)
#   source scripts/build_crc.csh
#
# What this does:
#   1. Loads the CRC modules AZURE2 needs (Python, Qt5, GSL, CMake, Qwt).
#   2. Installs pybind11/numpy into your own ~/.local if missing (needed to
#      build the pyazr `_azure2` extension).
#   3. Configures a Release build in ./build with the GUI, pyazr, MCMC,
#      ERYA (SRIM stopping powers), Qwt plotting and readline all enabled.
#   4. Builds with `make -jN`.
#   5. Copies the resulting executable to ~/bin/AZURE2 (creating ~/bin if
#      needed).
#
# Not needed on this cluster (already handled elsewhere):
#   - readline headers: system-wide (readline-devel), no module for it.
#   - Minuit2, the Coulomb-wavefunction library, and the SRIM stopping-power
#     utilities: vendored under external/, built in-tree automatically.
#   - GSL/Qt5/Qwt install *paths*: once the modules below are loaded, CMake
#     finds all three automatically via pkg-config / their CMake config
#     files -- no -DGSL_PATH, -DMINUIT_PATH or -DCMAKE_PREFIX_PATH needed.
#
# Safe to re-run any time (e.g. after `git pull`) to rebuild.

set BUILD_DIR = build

echo "Loading modules..."
module load python/3.12.13
module load qt/5.15.19
module load gsl/gcc/2.8
module load cmake/4.1.0
module load /groups/nsl/nslmodule/qwt/6.3.0/gcc/11.5.0

echo "Checking Python build deps (pybind11, numpy)..."
python3 -c "import pybind11" >& /dev/null
if ($status != 0) then
    echo "  installing pybind11 (--user)..."
    pip install --user pybind11
endif
python3 -c "import numpy" >& /dev/null
if ($status != 0) then
    echo "  installing numpy (--user)..."
    pip install --user numpy
endif
# scipy/mpmath aren't needed to *build* AZURE2, but pyazr-driven fit
# scripts import them transitively -- worth having up front.
python3 -c "import scipy, mpmath" >& /dev/null
if ($status != 0) then
    echo "  installing scipy, mpmath (--user, for pyazr scripts)..."
    pip install --user scipy mpmath
endif

set PYBIND11_DIR = `python3 -m pybind11 --cmakedir`
echo "pybind11 cmake dir: $PYBIND11_DIR"

echo "Configuring CMake ($BUILD_DIR)..."
cmake -S . -B $BUILD_DIR \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_GUI=ON \
    -DUSE_API=ON \
    -DUSE_QWT=ON \
    -DUSE_MCMC=ON \
    -DUSE_ERYA=ON \
    -DUSE_READLINE=ON \
    -DUSE_STAT=ON \
    -Dpybind11_DIR="$PYBIND11_DIR"
if ($status != 0) then
    echo "CMake configure failed -- see above."
    exit 1
endif

echo "Building..."
cmake --build $BUILD_DIR --parallel `nproc`
if ($status != 0) then
    echo "Build failed -- see above."
    exit 1
endif

echo "Installing to ~/bin/AZURE2..."
mkdir -p ~/bin
cp $BUILD_DIR/src/AZURE2 ~/bin/AZURE2

echo ""
echo "Done. Executable: ~/bin/AZURE2"
echo "If ~/bin isn't already on your PATH, add this to your .login:"
echo '    setenv PATH ${PATH}:${HOME}/bin'
echo ""
echo "pyazr lives in ./pyazr under this checkout. To 'import pyazr' from"
echo "elsewhere, either run scripts from this directory, or add this repo"
echo "to PYTHONPATH, e.g. in .login:"
echo '    setenv PYTHONPATH ${PYTHONPATH}:'"$cwd"
