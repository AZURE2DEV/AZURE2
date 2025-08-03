#!/bin/bash

# AZURE2 Linux Build Script
# This script builds AZURE2 for Linux systems

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}AZURE2 Linux Build Script${NC}"
echo "========================"

# Check if we're on Linux
if [[ "$OSTYPE" != "linux-gnu"* ]]; then
    echo -e "${RED}Error: This script must be run on Linux${NC}"
    exit 1
fi

# Create build directory
BUILD_DIR="build-linux"
rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR
cd $BUILD_DIR

echo -e "${GREEN}Configuring for Linux...${NC}"

# Auto-detect QWT installation
QWT_CMAKE_ARGS=""
if command -v pkg-config &> /dev/null && pkg-config --exists Qt5Qwt6 &> /dev/null; then
    echo -e "${YELLOW}Found QWT via pkg-config${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=ON"
elif [ -n "$CONDA_PREFIX" ] && [ -f "$CONDA_PREFIX/lib/libqwt.so" ]; then
    echo -e "${YELLOW}Found QWT in conda environment: $CONDA_PREFIX${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=ON -DCMAKE_PREFIX_PATH=$CONDA_PREFIX"
elif [ -f "/usr/local/lib/libqwt.so" ]; then
    echo -e "${YELLOW}Found QWT in /usr/local${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=ON -DCMAKE_PREFIX_PATH=/usr/local"
elif [ -f "/usr/lib/x86_64-linux-gnu/libqwt-qt5.so" ]; then
    echo -e "${YELLOW}Found QWT in system libraries${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=ON"
else
    echo -e "${YELLOW}QWT not found - building without plotting support${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=OFF"
fi

# Configure with CMake for Linux
cmake .. \
    -DBUILD_GUI=ON \
    -DUSE_API=ON \
    -DUSE_QWT=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    $QWT_CMAKE_ARGS

echo -e "${GREEN}Building AZURE2 for Linux...${NC}"
make -j$(nproc) VERBOSE=1

echo -e "${GREEN}Installing...${NC}"
make install

echo -e "${GREEN}Linux build completed successfully!${NC}"
echo "Executable: ${PWD}/src/AZURE2"
echo "GUI executable: ${PWD}/gui/AZURE2Setup"