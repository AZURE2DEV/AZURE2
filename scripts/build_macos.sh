#!/bin/bash

# AZURE2 macOS Application Bundle Build Script
# This script builds AZURE2 as a macOS application bundle (.app) and creates a DMG

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}AZURE2 macOS Bundle Build Script${NC}"
echo "================================"

# Check if we're on macOS
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo -e "${RED}Error: This script must be run on macOS${NC}"
    exit 1
fi

# Create build directory
BUILD_DIR="build-macos"
rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR
cd $BUILD_DIR

echo -e "${GREEN}Configuring for macOS bundle...${NC}"

# Auto-detect QWT installation
QWT_CMAKE_ARGS=""
if command -v pkg-config &> /dev/null && pkg-config --exists Qt5Qwt6 &> /dev/null; then
    echo -e "${YELLOW}Found QWT via pkg-config${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=ON"
elif [ -n "$CONDA_PREFIX" ] && [ -f "$CONDA_PREFIX/lib/libqwt.dylib" ]; then
    echo -e "${YELLOW}Found QWT in conda environment: $CONDA_PREFIX${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=ON -DCMAKE_PREFIX_PATH=$CONDA_PREFIX"
elif [ -f "/usr/local/lib/libqwt.dylib" ]; then
    echo -e "${YELLOW}Found QWT in /usr/local${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=ON -DCMAKE_PREFIX_PATH=/usr/local"
elif [ -f "/opt/homebrew/lib/libqwt.dylib" ]; then
    echo -e "${YELLOW}Found QWT in Homebrew${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=ON -DCMAKE_PREFIX_PATH=/opt/homebrew"
else
    echo -e "${YELLOW}QWT not found - building without plotting support${NC}"
    QWT_CMAKE_ARGS="-DUSE_QWT=OFF"
fi

# Configure with CMake for macOS bundle
cmake .. \
    -DBUILD_MACOS_BUNDLE=ON \
    -DBUILD_GUI=ON \
    -DUSE_API=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_OSX_DEPLOYMENT_TARGET=10.15 \
    -DCMAKE_INSTALL_PREFIX=/Applications \
    $QWT_CMAKE_ARGS

echo -e "${GREEN}Building AZURE2 for macOS...${NC}"
make -j$(sysctl -n hw.ncpu) VERBOSE=1

echo -e "${GREEN}Installing...${NC}"
make install

echo -e "${GREEN}Creating macOS application bundle...${NC}"

# Fix bundle structure if needed
BUNDLE_PATH="src/AZURE2.app"
if [ -d "$BUNDLE_PATH" ]; then
    echo "Bundle created at: $BUNDLE_PATH"
    
    # Copy any additional resources
    mkdir -p "$BUNDLE_PATH/Contents/Resources"
    
    # Copy SRIM2013.xml for ERYA support
    echo -e "${YELLOW}Copying SRIM2013.xml for ERYA support...${NC}"
    if [ -f "../erya/data/SRIM2013.xml" ]; then
        cp "../erya/data/SRIM2013.xml" "$BUNDLE_PATH/Contents/Resources/"
        echo "SRIM2013.xml copied to bundle Resources"
    elif [ -f "SRIM2013.xml" ]; then
        cp "SRIM2013.xml" "$BUNDLE_PATH/Contents/Resources/"
        echo "SRIM2013.xml copied to bundle Resources"
    else
        echo "Warning: SRIM2013.xml not found - ERYA features may not work"
    fi
    
    echo -e "${GREEN}Bundling dependencies for standalone distribution...${NC}"
    
    # Use macdeployqt if available (Qt's deployment tool)
    if command -v macdeployqt &> /dev/null; then
        echo -e "${YELLOW}Using macdeployqt to bundle Qt dependencies...${NC}"
        macdeployqt "$BUNDLE_PATH" -verbose=2 || echo -e "${YELLOW}macdeployqt failed, continuing with manual bundling${NC}"
    else
        echo -e "${YELLOW}macdeployqt not found, using manual dependency bundling${NC}"
    fi
    
    # Bundle additional dependencies manually
    echo -e "${YELLOW}Bundling additional dependencies...${NC}"
    
    FRAMEWORKS_DIR="$BUNDLE_PATH/Contents/Frameworks"
    LIBS_DIR="$BUNDLE_PATH/Contents/Libraries"
    mkdir -p "$FRAMEWORKS_DIR" "$LIBS_DIR"
    
    # Function to copy and fix library paths
    bundle_library() {
        local lib_path="$1"
        local lib_name=$(basename "$lib_path")
        
        if [ -f "$lib_path" ] && [ ! -f "$LIBS_DIR/$lib_name" ]; then
            echo "Bundling: $lib_name"
            cp "$lib_path" "$LIBS_DIR/"
            
            # Fix the library's internal ID and dependency paths
            install_name_tool -id "@loader_path/../Libraries/$lib_name" "$LIBS_DIR/$lib_name" 2>/dev/null || true
            
            # Update the main executable to use the bundled library
            install_name_tool -change "$lib_path" "@loader_path/../Libraries/$lib_name" "$BUNDLE_PATH/Contents/MacOS/AZURE2" 2>/dev/null || true
            install_name_tool -change "/usr/local/lib/$lib_name" "@loader_path/../Libraries/$lib_name" "$BUNDLE_PATH/Contents/MacOS/AZURE2" 2>/dev/null || true
            install_name_tool -change "/opt/homebrew/lib/$lib_name" "@loader_path/../Libraries/$lib_name" "$BUNDLE_PATH/Contents/MacOS/AZURE2" 2>/dev/null || true
            
            # Also fix other bundled libraries that might depend on this one
            for bundled_lib in "$LIBS_DIR"/*.dylib; do
                if [ -f "$bundled_lib" ]; then
                    install_name_tool -change "$lib_path" "@loader_path/$lib_name" "$bundled_lib" 2>/dev/null || true
                    install_name_tool -change "/usr/local/lib/$lib_name" "@loader_path/$lib_name" "$bundled_lib" 2>/dev/null || true
                    install_name_tool -change "/opt/homebrew/lib/$lib_name" "@loader_path/$lib_name" "$bundled_lib" 2>/dev/null || true
                fi
            done
        fi
    }
    
    # Get all dynamic library dependencies
    echo -e "${YELLOW}Analyzing dependencies...${NC}"
    DEPENDENCIES=$(otool -L "$BUNDLE_PATH/Contents/MacOS/AZURE2" | grep -E "\t/" | grep -v "/System/Library\|/usr/lib\|@" | awk '{print $1}')
    
    # Bundle each dependency
    for dep in $DEPENDENCIES; do
        bundle_library "$dep"
        
        # Also bundle dependencies of dependencies (one level deep)
        if [ -f "$dep" ]; then
            SUB_DEPS=$(otool -L "$dep" 2>/dev/null | grep -E "\t/" | grep -v "/System/Library\|/usr/lib\|@" | awk '{print $1}')
            for sub_dep in $SUB_DEPS; do
                bundle_library "$sub_dep"
            done
        fi
    done
    
    # Special handling for common library locations
    for lib_dir in "/usr/local/lib" "/opt/homebrew/lib" "$CONDA_PREFIX/lib"; do
        if [ -d "$lib_dir" ]; then
            # Bundle GSL libraries
            for gsl_lib in "$lib_dir"/libgsl*.dylib; do
                [ -f "$gsl_lib" ] && bundle_library "$gsl_lib"
            done
            
            # Bundle QWT library
            for qwt_lib in "$lib_dir"/libqwt*.dylib; do
                [ -f "$qwt_lib" ] && bundle_library "$qwt_lib"
            done
            
            # Bundle other common scientific libraries
            for lib in "$lib_dir"/lib{blas,lapack,openblas,readline,ncurses}*.dylib; do
                [ -f "$lib" ] && bundle_library "$lib"
            done
        fi
    done
    
    # Fix any remaining absolute paths in bundled libraries
    echo -e "${YELLOW}Fixing library paths...${NC}"
    for lib in "$LIBS_DIR"/*.dylib; do
        if [ -f "$lib" ]; then
            # Get all dependencies of this library
            LIB_DEPS=$(otool -L "$lib" 2>/dev/null | grep -E "\t/" | grep -v "/System/Library\|/usr/lib\|@" | awk '{print $1}')
            for lib_dep in $LIB_DEPS; do
                lib_dep_name=$(basename "$lib_dep")
                if [ -f "$LIBS_DIR/$lib_dep_name" ]; then
                    install_name_tool -change "$lib_dep" "@loader_path/$lib_dep_name" "$lib" 2>/dev/null || true
                fi
            done
        fi
    done
    
    # Verify the bundle
    echo -e "${YELLOW}Verifying bundle dependencies...${NC}"
    if otool -L "$BUNDLE_PATH/Contents/MacOS/AZURE2" | grep -q "/usr/local\|/opt/homebrew" | grep -v "@"; then
        echo -e "${YELLOW}Warning: Some external dependencies may still be present${NC}"
        otool -L "$BUNDLE_PATH/Contents/MacOS/AZURE2" | grep "/usr/local\|/opt/homebrew" | head -5
    else
        echo -e "${GREEN}All external dependencies appear to be bundled${NC}"
    fi
    
    echo -e "${GREEN}Dependency bundling completed${NC}"
    echo "Bundled libraries: $(ls -1 "$LIBS_DIR"/*.dylib 2>/dev/null | wc -l)"
    
    # Sign the bundle (if developer certificates are available)
    if command -v codesign &> /dev/null; then
        echo -e "${YELLOW}Attempting to sign the application bundle...${NC}"
        codesign --deep --force --verify --verbose --sign - "$BUNDLE_PATH" || echo -e "${YELLOW}Code signing failed (no certificates available)${NC}"
    fi
    
    echo -e "${GREEN}Creating DMG package...${NC}"
    make package
    
    # Create enhanced distributable DMG with documentation
    if [ -f "AZURE2-*-macOS.dmg" ]; then
        echo -e "${YELLOW}Creating enhanced DMG with documentation...${NC}"
        
        # Create temporary directory for DMG contents
        DMG_TEMP=$(mktemp -d)
        mkdir -p "$DMG_TEMP/AZURE2"
        
        # Copy the application bundle
        cp -R "$BUNDLE_PATH" "$DMG_TEMP/AZURE2/"
        
        # Copy documentation if available
        if [ -f "../README.md" ]; then
            cp "../README.md" "$DMG_TEMP/AZURE2/"
        fi
        
        if [ -f "../LICENSE" ]; then
            cp "../LICENSE" "$DMG_TEMP/AZURE2/"
        fi
        
        # Create README for distribution
        cat > "$DMG_TEMP/AZURE2/README_DISTRIBUTION.txt" << EOF

AZURE2 for macOS - Standalone Distribution
==========================================

This is a standalone distribution of AZURE2 that includes all necessary dependencies.

Installation:
1. Drag AZURE2.app to your Applications folder (optional)
2. Double-click AZURE2.app to run

System Requirements:
- macOS 10.15 (Catalina) or later
- 64-bit Intel or Apple Silicon Mac

All required libraries have been bundled within the application.
No additional software installation is required.

For documentation and support, visit: https://azure.nd.edu/

Build Information:
- Built on: $(date)
- Bundled libraries: $(ls -1 "$LIBS_DIR"/*.dylib 2>/dev/null | wc -l)
- Architecture: $(uname -m)
EOF
        
        # Create a symbolic link to Applications for easy installation
        ln -sf /Applications "$DMG_TEMP/Applications"
        
        # Create the final DMG
        FINAL_DMG_NAME="AZURE2-$(date +%Y%m%d)-macOS-Standalone.dmg"
        hdiutil create -srcfolder "$DMG_TEMP" -volname "AZURE2" -fs HFS+ -fsargs "-c c=64,a=16,e=16" -format UDBZ "$FINAL_DMG_NAME"
        
        # Cleanup
        rm -rf "$DMG_TEMP"
        
        echo -e "${GREEN}Enhanced DMG created: $FINAL_DMG_NAME${NC}"
    fi
    
    echo -e "${GREEN}macOS build completed successfully!${NC}"
    echo "Application bundle: ${PWD}/${BUNDLE_PATH}"
    echo "Standard DMG: ${PWD}/AZURE2-*-macOS.dmg"
    echo "Standalone DMG: ${PWD}/AZURE2-*-macOS-Standalone.dmg"
    
    # Show final bundle information
    echo -e "${BLUE}Final bundle information:${NC}"
    echo "Bundle size: $(du -sh "$BUNDLE_PATH" | cut -f1)"
    echo "Bundled dependencies: $(ls -1 "$LIBS_DIR"/*.dylib 2>/dev/null | wc -l || echo 0)"
    
    # Verify the bundle can run (basic check)
    if [ -x "$BUNDLE_PATH/Contents/MacOS/AZURE2" ]; then
        echo -e "${GREEN}Bundle executable is present and executable${NC}"
    else
        echo -e "${RED}Warning: Bundle executable may have issues${NC}"
    fi
else
    echo -e "${RED}Error: Application bundle was not created${NC}"
    exit 1
fi