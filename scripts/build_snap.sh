#!/bin/bash

# AZURE2 Snap Package Build Script
# This script builds AZURE2 as a snap package for Linux distribution

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}AZURE2 Snap Package Build Script${NC}"
echo "================================="

# Check if snapcraft is installed
if ! command -v snapcraft &> /dev/null; then
    echo -e "${RED}Error: snapcraft is not installed${NC}"
    echo "To install snapcraft:"
    echo "  Ubuntu/Debian: sudo snap install snapcraft --classic"
    echo "  Fedora: sudo dnf install snapcraft"
    echo "  Arch: yay -S snapcraft"
    exit 1
fi

# Check if running on macOS
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo -e "${YELLOW}Warning: Running on macOS${NC}"
    echo "Snapcraft on macOS requires Linux VM and may have file permission issues."
    echo "Alternative approaches:"
    echo "1. Use Docker to build snap in Linux container"
    echo "2. Use GitHub Actions or CI/CD for snap building"
    echo "3. Build on actual Linux system"
    echo ""
    echo "Attempting Docker-based build instead..."
    
    # Check if Docker is available
    if command -v docker &> /dev/null; then
        echo -e "${BLUE}Using Docker to build snap package...${NC}"
        
        # Use dedicated Docker container for snap building
        if [ -f "docker/Dockerfile.snap" ]; then
            echo "Building snap using dedicated Docker container..."
            docker build -f docker/Dockerfile.snap -t azure2-snap-builder --platform linux/amd64 .
            docker run --rm -v "$(pwd)":/output -w /output azure2-snap-builder bash -c "
                cp /build/azure2_*.snap . 2>/dev/null || echo 'No snap file found'
                ls -la azure2_*.snap 2>/dev/null || echo 'Build may have failed'
            "
        else
            # Fallback to inline Docker build
            echo "Using inline Docker build..."
            docker run --rm -v "$(pwd)":/build -w /build ubuntu:22.04 bash -c "
                apt-get update && 
                apt-get install -y snapcraft build-essential cmake pkg-config libgsl-dev qtbase5-dev libqwt-qt5-dev && 
                chmod 755 snap && chmod 644 snap/snapcraft.yaml &&
                snapcraft --verbose
            "
        fi
        exit $?
    else
        echo -e "${RED}Docker not available. Please install Docker or use Linux system.${NC}"
        exit 1
    fi
fi

# Check if snap directory exists
if [ ! -d "snap" ]; then
    echo -e "${RED}Error: snap directory not found${NC}"
    echo "Please run this script from the AZURE2 root directory"
    exit 1
fi

# Check if snapcraft.yaml exists
if [ ! -f "snap/snapcraft.yaml" ]; then
    echo -e "${RED}Error: snap/snapcraft.yaml not found${NC}"
    echo "The snap configuration file is missing"
    exit 1
fi

echo -e "${GREEN}Building AZURE2 snap package...${NC}"

# Clean any previous builds
echo -e "${YELLOW}Cleaning previous builds...${NC}"
snapcraft clean || true

# Build the snap
echo -e "${YELLOW}Starting snap build process...${NC}"
echo "This may take several minutes..."

# Build with detailed output
if snapcraft --verbose --debug; then
    echo -e "${GREEN}Snap build completed successfully!${NC}"
    
    # Find the generated snap file
    SNAP_FILE=$(find . -name "azure2_*.snap" -type f | head -1)
    
    if [ -n "$SNAP_FILE" ]; then
        echo -e "${GREEN}Generated snap package: $SNAP_FILE${NC}"
        
        # Show file information
        echo -e "${BLUE}Package information:${NC}"
        ls -lh "$SNAP_FILE"
        
        # Show package details
        echo -e "${BLUE}Package details:${NC}"
        snap info "$SNAP_FILE" 2>/dev/null || echo "Run 'snap info $SNAP_FILE' to see package details"
        
        echo -e "${GREEN}Installation instructions:${NC}"
        echo "To install locally:"
        echo "  sudo snap install $SNAP_FILE --dangerous"
        echo ""
        echo "To install for development/testing:"
        echo "  sudo snap install $SNAP_FILE --dangerous --devmode"
        echo ""
        echo "After installation, run with:"
        echo "  azure2        # Command line interface"
        echo "  azure2.azure2-gui  # Graphical interface"
        echo ""
        echo "To publish to the Snap Store:"
        echo "  snapcraft upload $SNAP_FILE"
        echo "  snapcraft release azure2 <revision> stable"
        
    else
        echo -e "${RED}Error: Snap file was not created${NC}"
        exit 1
    fi
    
else
    echo -e "${RED}Snap build failed!${NC}"
    echo "Common issues and solutions:"
    echo "1. Missing dependencies: Install required development packages"
    echo "2. Network issues: Check internet connection for downloading packages"
    echo "3. Disk space: Ensure sufficient disk space for build"
    echo "4. Snapcraft version: Try updating snapcraft"
    echo ""
    echo "For detailed logs, check the output above or run:"
    echo "  snapcraft --debug"
    exit 1
fi

echo -e "${GREEN}Snap packaging completed!${NC}"