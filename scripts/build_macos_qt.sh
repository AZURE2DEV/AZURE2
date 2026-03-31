#!/bin/bash
# AZURE2 macOS Application Bundle Build Script
# Automatically detects Qt6 + QWT, builds the .app, bundles dependencies, and creates DMG.

set -e

# --- LLVM / Clang setup ---
export CC=/opt/homebrew/opt/llvm/bin/clang
export CXX=/opt/homebrew/opt/llvm/bin/clang++
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"

# --- Qt + QWT setup (Homebrew) ---
export PATH="$(brew --prefix qt)/bin:$PATH"
export CMAKE_PREFIX_PATH="$(brew --prefix qt):$(brew --prefix qwt)"
export PKG_CONFIG_PATH="$(brew --prefix qwt)/lib/pkgconfig:$PKG_CONFIG_PATH"

# --- Colors for terminal output ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${GREEN}AZURE2 macOS Bundle Build Script${NC}"
echo "================================"

# --- Check macOS ---
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo -e "${RED}Error: This script must be run on macOS${NC}"
    exit 1
fi

# --- Create clean build directory ---
BUILD_DIR="build-macos"
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

echo -e "${GREEN}Configuring CMake for macOS bundle...${NC}"

# --- Run CMake ---
cmake .. \
  -DBUILD_MACOS_BUNDLE=ON \
  -DBUILD_GUI=ON \
  -DUSE_API=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.15 \
  -DCMAKE_INSTALL_PREFIX=/Applications \
  -DUSE_QWT=ON \
  -DQWT_INCLUDE_DIR=$(brew --prefix qwt)/include \
  -DQWT_LIBRARY=$(brew --prefix qwt)/lib/qwt.framework/qwt

echo -e "${GREEN}Building AZURE2...${NC}"
make -j$(sysctl -n hw.ncpu) VERBOSE=1

echo -e "${GREEN}Installing...${NC}"
make install

# --- Bundle app ---
BUNDLE_PATH="src/AZURE2.app"
if [ ! -d "$BUNDLE_PATH" ]; then
    echo -e "${RED}Error: Application bundle was not created${NC}"
    exit 1
fi

echo "Bundle created at: $BUNDLE_PATH"

# --- Copy additional resources ---
mkdir -p "$BUNDLE_PATH/Contents/Resources"
if [ -f "../erya/data/SRIM2013.xml" ]; then
    cp "../erya/data/SRIM2013.xml" "$BUNDLE_PATH/Contents/Resources/"
    echo "SRIM2013.xml copied to bundle Resources"
elif [ -f "SRIM2013.xml" ]; then
    cp "SRIM2013.xml" "$BUNDLE_PATH/Contents/Resources/"
    echo "SRIM2013.xml copied to bundle Resources"
else
    echo "Warning: SRIM2013.xml not found - ERYA features may not work"
fi

# --- macdeployqt (if available) ---
if command -v macdeployqt &> /dev/null; then
    echo -e "${YELLOW}Bundling Qt frameworks using macdeployqt...${NC}"
    macdeployqt "$BUNDLE_PATH" -verbose=2 || echo -e "${YELLOW}macdeployqt failed, continuing${NC}"
else
    echo -e "${YELLOW}macdeployqt not found; you may need to manually bundle Qt${NC}"
fi

# --- Optional: manual dependency bundling (GSL, OpenMP, QWT, etc.) ---
LIBS_DIR="$BUNDLE_PATH/Contents/Libraries"
mkdir -p "$LIBS_DIR"

bundle_library() {
    local lib="$1"
    local name=$(basename "$lib")
    if [ -f "$lib" ] && [ ! -f "$LIBS_DIR/$name" ]; then
        echo "Bundling $name"
        cp "$lib" "$LIBS_DIR/"
        install_name_tool -id "@loader_path/../Libraries/$name" "$LIBS_DIR/$name" 2>/dev/null || true
        install_name_tool -change "$lib" "@loader_path/../Libraries/$name" "$BUNDLE_PATH/Contents/MacOS/AZURE2" 2>/dev/null || true
    fi
}

# Bundle GSL, QWT, etc.
for lib_dir in "/opt/homebrew/lib" "$CONDA_PREFIX/lib"; do
    [ -d "]()
