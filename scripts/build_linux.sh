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

# Configure with CMake for Linux
cmake .. \
    -DUSE_API=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr/local

echo -e "${GREEN}Building AZURE2 for Linux...${NC}"
make -j$(nproc) VERBOSE=1

echo -e "${GREEN}Installing...${NC}"
make install

echo -e "${GREEN}Linux build completed successfully!${NC}"
echo "Executable: ${PWD}/src/AZURE2"