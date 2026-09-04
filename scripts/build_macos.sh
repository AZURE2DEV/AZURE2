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
    -DCMAKE_INSTALL_PREFIX="$PWD/stage" \
    $QWT_CMAKE_ARGS

echo -e "${GREEN}Building AZURE2 for macOS...${NC}"
make -j$(sysctl -n hw.ncpu) VERBOSE=1

# Into a staging directory under the build tree, never into /Applications
# itself.  This script's deliverable is the .dmg; the user drags the app across
# from it.  Installing straight to the prefix also scattered the project's
# headers, static libraries and cmake/pkg-config files through
# /Applications/include, /Applications/lib and /Applications/share, which is
# not what anyone means by "install an app".
echo -e "${GREEN}Staging install tree...${NC}"
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
    
    # Sweep the whole bundle to a fixpoint.  The pass above follows
    # dependencies one level deep, which is not enough -- readline pulls
    # ncurses and ncurses pulls libtinfo, so a second-level library was left
    # pointing at the build machine and the bundle failed on any Mac that did
    # not happen to have it.  Every Mach-O file is examined, the Qt plugins
    # included, since a plugin can need a library the executable never
    # references.
    echo -e "${YELLOW}Sweeping transitive dependencies...${NC}"

    bundle_machos() {
        find "$BUNDLE_PATH" -type f \( -name '*.dylib' -o -name '*.so' -o -perm +111 \) 2>/dev/null \
            | while read -r f; do file "$f" 2>/dev/null | grep -q 'Mach-O' && echo "$f"; done
    }
    # Where Contents/Libraries sits relative to a given file, so a plugin two
    # directories down and the executable one down each get a correct path.
    rel_to_libs() {
        python3 -c 'import os,sys; print(os.path.relpath(sys.argv[1], os.path.dirname(sys.argv[2])))' \
            "$LIBS_DIR" "$1"
    }

    for pass in 1 2 3 4 5 6; do
        added=0
        for f in $(bundle_machos); do
            install_name_tool -add_rpath "@loader_path/$(rel_to_libs "$f")" "$f" 2>/dev/null || true
            for dep in $(otool -L "$f" 2>/dev/null | tail -n +2 | awk '{print $1}'); do
                case "$dep" in @*|/usr/lib/*|/System/*) continue ;; esac
                name=$(basename "$dep")
                if [ ! -f "$LIBS_DIR/$name" ] && [ -f "$dep" ]; then
                    cp -L "$dep" "$LIBS_DIR/$name" 2>/dev/null || continue
                    chmod u+w "$LIBS_DIR/$name"
                    install_name_tool -id "@rpath/$name" "$LIBS_DIR/$name" 2>/dev/null || true
                    echo "  + $name"
                    added=1
                fi
                [ -f "$LIBS_DIR/$name" ] && install_name_tool -change "$dep" "@rpath/$name" "$f" 2>/dev/null || true
            done
        done
        if [ "$added" -eq 0 ]; then echo "  settled after pass $pass"; break; fi
    done

    # Verify for real: anything still naming an absolute path outside the
    # bundle is a library this Mac happens to have and the user's may not.
    echo -e "${YELLOW}Verifying bundle dependencies...${NC}"
    # Collected through a file rather than $(...): bash 3.2, which is what
    # /bin/bash still is on macOS, mis-parses a case pattern's ")" inside a
    # command substitution.
    LEAKFILE=$(mktemp)
    for f in $(bundle_machos); do
        for d in $(otool -L "$f" 2>/dev/null | tail -n +2 | awk '{print $1}'); do
            case "$d" in
                @*|/usr/lib/*|/System/*) continue ;;
            esac
            echo "$(basename "$f") -> $d" >> "$LEAKFILE"
        done
    done
    sort -u "$LEAKFILE" -o "$LEAKFILE"
    if [ -s "$LEAKFILE" ]; then
        echo -e "${RED}The bundle still references libraries outside itself:${NC}"
        cat "$LEAKFILE"; rm -f "$LEAKFILE"
        echo -e "${RED}It would fail on a Mac without them.${NC}"
        exit 1
    fi
    rm -f "$LEAKFILE"
    if [ ! -f "$BUNDLE_PATH/Contents/PlugIns/platforms/libqcocoa.dylib" ]; then
        echo -e "${RED}No Contents/PlugIns/platforms/libqcocoa.dylib -- Qt cannot start${NC}"
        echo -e "${RED}without its platform plugin and the app dies before main().${NC}"
        exit 1
    fi
    echo -e "${GREEN}Self-contained: no external references, platform plugin present${NC}"

    echo -e "${GREEN}Dependency bundling completed${NC}"
    echo "Bundled libraries: $(ls -1 "$LIBS_DIR"/*.dylib 2>/dev/null | wc -l)"
    
    # Sign the bundle (if developer certificates are available)
    if command -v codesign &> /dev/null; then
        echo -e "${YELLOW}Attempting to sign the application bundle...${NC}"
        codesign --deep --force --verify --verbose --sign - "$BUNDLE_PATH" || echo -e "${YELLOW}Code signing failed (no certificates available)${NC}"
    fi
    
    # The disk image is built here rather than with "make package".  CPack's
    # DragNDrop generator packages the whole install tree, so its .dmg carried
    # include/, lib/ and share/ -- the project's headers and static libraries --
    # beside the application, which is confusing to hand to a student.  This one
    # holds the app, an Applications symlink to drag it onto, and a note.
    echo -e "${GREEN}Creating DMG package...${NC}"
    if true; then
        DMG_TEMP=$(mktemp -d)

        # At the top level, beside the Applications symlink: the drag-and-drop
        # gesture the window is asking for only reads if the two sit together.
        cp -R "$BUNDLE_PATH" "$DMG_TEMP/"
        
        # Copy documentation if available
        if [ -f "../LICENSE" ]; then
            cp "../LICENSE" "$DMG_TEMP/LICENSE.txt"
        fi

        # Create README for distribution
        cat > "$DMG_TEMP/READ ME FIRST.txt" << EOF

AZURE2 for macOS - Standalone Distribution
==========================================

This is a standalone distribution of AZURE2 that includes all necessary dependencies.

Installation:
1. Drag AZURE2.app onto the Applications folder shown beside it.
2. The first time you open it, macOS will say the developer cannot be
   verified: this build is signed ad-hoc, not notarized with an Apple
   Developer ID. Right-click (or control-click) AZURE2.app, choose "Open",
   and confirm. Only the first launch needs this.

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
        FINAL_DMG_NAME="AZURE2-$(date +%Y%m%d)-macOS-$(uname -m).dmg"
        hdiutil create -srcfolder "$DMG_TEMP" -volname "AZURE2" -fs HFS+ -fsargs "-c c=64,a=16,e=16" -format UDBZ "$FINAL_DMG_NAME"
        
        # Cleanup
        rm -rf "$DMG_TEMP"
        
        echo -e "${GREEN}Enhanced DMG created: $FINAL_DMG_NAME${NC}"
    fi
    
    echo -e "${GREEN}macOS build completed successfully!${NC}"
    echo "Application bundle: ${PWD}/${BUNDLE_PATH}"
    echo "Disk image:        ${PWD}/${FINAL_DMG_NAME}"
    
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