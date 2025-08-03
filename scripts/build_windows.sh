#!/bin/bash

# AZURE2 Windows Docker Build Script
# This script builds AZURE2 for Windows using Docker with MinGW cross-compilation

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}AZURE2 Windows Docker Build Script${NC}"
echo "=================================="

# Check if Docker is installed and running
if ! command -v docker &> /dev/null; then
    echo -e "${RED}Error: Docker is not installed!${NC}"
    echo "Please install Docker:"
    echo "  macOS: brew install --cask docker"
    echo "  Linux: Follow Docker installation guide for your distribution"
    exit 1
fi

if ! docker info &> /dev/null; then
    echo -e "${RED}Error: Docker is not running!${NC}"
    echo "Please start Docker and try again."
    exit 1
fi

# Setup output directory
OUTPUT_DIR="./build-windows"
# If existing, remove the old output directory
if [ -d "$OUTPUT_DIR" ]; then
    echo -e "${YELLOW}Removing old output directory: $OUTPUT_DIR${NC}"
    rm -rf "$OUTPUT_DIR"
fi

# Check if Docker image needs to be rebuilt
SKIP_BUILD=false
if [ "$SKIP_BUILD" = false ]; then
    # Build image (targeting builder stage to get build environment)
    echo -e "${GREEN}Building Windows compilation Docker image...${NC}"
    echo -e "${YELLOW}Building with GUI and QWT support for Windows${NC}"
    docker build -f docker/Dockerfile.windows --target builder -t azure2-windows-builder --label dockerfile_hash="$DOCKERFILE_HASH" . 2>&1 | tee build.log

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo -e "${RED}Docker build failed! Check build.log for details.${NC}"
        echo -e "${YELLOW}Common issues and solutions:${NC}"
        echo "1. Missing packages: Check if all mingw64-* packages are available"
        echo "2. CMake configuration: Verify CMakeLists.txt cross-compilation settings"
        echo "3. Library linking: Check Qt5 and QWT library paths"
        echo "4. Compilation errors: Review C++ source compatibility with MinGW"
        exit 1
    fi
else
    echo -e "${GREEN}Using existing Docker image${NC}"
fi

# Extract the Windows executable
echo -e "${GREEN}Extracting Windows executable...${NC}"
# Create a new output directory
echo -e "${YELLOW}Creating output directory: $OUTPUT_DIR${NC}"
mkdir -p "$OUTPUT_DIR"

# Copy the executable and DLLs from the Docker image
echo -e "${YELLOW}Copying executable and DLLs from Docker container...${NC}"

# Create container from the built image (defaults to final stage)
CONTAINER_ID=$(docker create azure2-windows-builder)

# Start container to execute commands
docker start "$CONTAINER_ID" >/dev/null 2>&1

# Try to find the executable in different possible locations
EXE_FOUND=false
for exe_path in "/workspace/build-windows/src/AZURE2.exe" "/workspace/build-windows/AZURE2.exe" "/usr/local/bin/AZURE2.exe"; do
    if docker exec "$CONTAINER_ID" test -f "$exe_path" 2>/dev/null; then
        echo -e "${GREEN}Found executable at: $exe_path${NC}"
        docker cp "$CONTAINER_ID:$exe_path" "$OUTPUT_DIR/" 2>/dev/null && EXE_FOUND=true
        break
    fi
done

if [ "$EXE_FOUND" = false ]; then
    echo -e "${RED}Error: Could not find AZURE2.exe in any expected location${NC}"
    echo "All files in build directory:"
    docker exec "$CONTAINER_ID" find /workspace/build-windows -type f 2>/dev/null | head -20
    docker stop "$CONTAINER_ID" >/dev/null 2>&1
    docker rm "$CONTAINER_ID"
    exit 1
fi

# Copy DLL files from multiple possible locations
echo -e "${YELLOW}Copying DLL files...${NC}"

# List specific DLL files that should be copied
EXPECTED_DLLS="libgsl-0.dll libgslcblas-0.dll libreadline8.dll libhistory8.dll Qt5Core.dll Qt5Gui.dll Qt5Widgets.dll Qt5PrintSupport.dll Qt5OpenGL.dll Qt5Svg.dll qwt-qt5.dll libgcc_s_seh-1.dll libstdc++-6.dll libwinpthread-1.dll iconv.dll libiconv-2.dll libpcre2-16-0.dll libpcre2-8-0.dll zlib1.dll libharfbuzz-0.dll libpng16-16.dll libfreetype-6.dll libglib-2.0-0.dll libintl-8.dll libgraphite2.dll libbz2-1.dll libtermcap-0.dll"

for dll in $EXPECTED_DLLS; do
    DLL_FOUND=false
    # Check multiple possible locations for DLLs
    for dll_path in "/workspace/build-windows/src/$dll" "/workspace/build-windows/$dll" "/usr/local/bin/$dll"; do
        if docker exec "$CONTAINER_ID" test -f "$dll_path" 2>/dev/null; then
            echo "Copying DLL: $dll from $dll_path"
            docker cp "$CONTAINER_ID:$dll_path" "$OUTPUT_DIR/" 2>/dev/null && DLL_FOUND=true
            break
        fi
    done
    
    if [ "$DLL_FOUND" = false ]; then
        echo "DLL not found in any location: $dll"
    fi
done

# Copy platforms directory from multiple possible locations
PLATFORMS_FOUND=false
for platforms_path in "/workspace/build-windows/src/platforms" "/workspace/build-windows/platforms" "/usr/local/bin/platforms"; do
    if docker exec "$CONTAINER_ID" test -d "$platforms_path" 2>/dev/null; then
        echo -e "${YELLOW}Copying Qt5 platforms directory from: $platforms_path${NC}"
        mkdir -p "$OUTPUT_DIR/platforms"
        docker cp "$CONTAINER_ID:$platforms_path/." "$OUTPUT_DIR/platforms/" 2>/dev/null && PLATFORMS_FOUND=true
        break
    fi
done

if [ "$PLATFORMS_FOUND" = false ]; then
    echo "Platforms directory not found in any expected location"
fi

# Cleanup
echo -e "${YELLOW}Cleaning up Docker image...${NC}"

# Stop and remove container
docker stop "$CONTAINER_ID" >/dev/null 2>&1
docker rm "$CONTAINER_ID"