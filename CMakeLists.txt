# Range form: needs at least 3.16, and opts in to the current behaviour of
# anything up to 4.0.  Requiring 4.0 outright kept the project from configuring
# anywhere but a very recent toolchain -- Ubuntu 22.04 ships CMake 3.22 and
# 24.04 ships 3.28, so distro builds could never get past line 1.
cmake_minimum_required(VERSION 3.16...4.0)
project(AZURE2)

# Default to an optimized build.  With no CMAKE_BUILD_TYPE, CMake passes no
# optimization flags at all, which for a numerical code means an interpreter-
# speed binary that looks like a normal one.  Multi-config generators (Visual
# Studio, Xcode) pick the configuration at build time, so leave those alone.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo MinSizeRel)
endif()

# Version info
set(AZURE2_VERSION_MAJOR 1)
set(AZURE2_VERSION_MINOR 1)
set(AZURE2_VERSION_PATCH 0)
set(AZURE2_VERSION "${AZURE2_VERSION_MAJOR}.${AZURE2_VERSION_MINOR}.${AZURE2_VERSION_PATCH}")

# M_PI and friends are not ISO C++, so MinGW's headers hide them whenever
# __STRICT_ANSI__ is set -- which -std=c++17 (as opposed to -std=gnu++17) does.
# _USE_MATH_DEFINES is the MSVC-originated convention for asking for them back,
# honoured by both MSVC and mingw-w64. It has to be visible before any header
# reaches <cmath>, hence a command-line definition rather than a #define in a
# source file. Also covers CROSS_COMPILE_WINDOWS, which sets CMAKE_SYSTEM_NAME.
if(WIN32)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()

# Options
option(CROSS_COMPILE_WINDOWS "Cross-compile for Windows" OFF)
option(BUILD_MACOS_BUNDLE "Build macOS application bundle" OFF)
option(BUILD_LIBRARY "Build AZURE2 as a library" OFF)
option(BUILD_GUI "Build AZURE2 with graphical setup program" ON)
option(USE_STAT "Use stat() function (not for Windows)" ON)
option(USE_QWT "Use Qwt to build Plot tab." ON)
option(USE_API "Use API version." ON)
if(USE_API)
    # The _azure2 pybind11 module (api/src/CMakeLists.txt) is a shared object
    # that statically links AZUREAPI plus its dependencies (AZURE2Core,
    # Minuit2, coul). All of those must be position-independent for the
    # link to succeed, so default every target to PIC before their
    # subdirectories are added.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()
option(USE_READLINE "Use Readline library for CLI input." ON)
option(USE_ERYA "Use ERYA SRIM stopping power utilities." ON)
option(USE_MCMC "Use MCMC Bayesian sampling OFF" ON)
option(USE_GCOV "Create a GCov-enabled build." OFF)
option(CODE_COVERAGE "Enable coverage reporting" OFF)

# ---------------------------------------------------------------- macOS/Homebrew
# Homebrew needs help that a Linux build never does, and without it a plain
# "cmake -S . -B build" fails on macOS for two unrelated reasons:
#
#   * Apple's clang ships no OpenMP runtime at all.  The compiler needs
#     "-Xpreprocessor -fopenmp", the runtime comes from the separate libomp
#     formula, and FindOpenMP cannot infer either for AppleClang -- so
#     find_package(OpenMP REQUIRED) below fails outright.
#   * qt@5 and qwt-qt5 are keg-only: Homebrew deliberately leaves them out of
#     the default prefix, so find_package(Qt5) does not see them.
#
# Everything here is derived by asking brew, and every value is overridable
# from the command line, so an explicit -D or a conda/MacPorts toolchain (which
# has none of these problems, having its own compiler and Qt in one prefix)
# still wins.
if(APPLE AND NOT CROSS_COMPILE_WINDOWS)
    find_program(BREW_EXECUTABLE brew)
    if(BREW_EXECUTABLE)
        # Keg-only or not, adding the prefix is harmless and makes the build
        # independent of what the user has linked into /opt/homebrew.
        foreach(_formula qt@5 qwt-qt5 qwt readline gsl libomp)
            execute_process(COMMAND ${BREW_EXECUTABLE} --prefix ${_formula}
                            OUTPUT_VARIABLE _prefix
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            ERROR_QUIET
                            RESULT_VARIABLE _rc)
            if(_rc EQUAL 0 AND IS_DIRECTORY "${_prefix}")
                list(APPEND CMAKE_PREFIX_PATH "${_prefix}")
                if(_formula STREQUAL "libomp")
                    set(_libomp_prefix "${_prefix}")
                endif()
                message(STATUS "Homebrew ${_formula}: ${_prefix}")
            endif()
        endforeach()

        # The libomp prefix is only recorded here.  It is applied further down,
        # and then only if the compiler cannot do OpenMP by itself -- a conda or
        # MacPorts toolchain brings its own, and must not be overridden.
    else()
        message(STATUS "Homebrew not found; assuming a self-contained toolchain "
                       "(conda/MacPorts) that provides OpenMP and Qt5 itself.")
    endif()
endif()

# Cross-compilation settings
if(CROSS_COMPILE_WINDOWS)
    set(CMAKE_SYSTEM_NAME Windows)
    set(CMAKE_C_COMPILER /usr/bin/x86_64-w64-mingw32-gcc)
    set(CMAKE_CXX_COMPILER /usr/bin/x86_64-w64-mingw32-g++)
    set(CMAKE_RC_COMPILER /usr/bin/x86_64-w64-mingw32-windres)
    set(CMAKE_FIND_ROOT_PATH /usr/x86_64-w64-mingw32)
    set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
    set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
    set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

    # Configure proper static linking for MinGW runtime libraries
    # Note: Full static linking with Qt5 is not possible, so we use selective static linking
    set(CMAKE_C_FLAGS "-static-libgcc")
    set(CMAKE_CXX_FLAGS "-static-libgcc -static-libstdc++")
    
    # Use advanced static linking approach for winpthread while allowing dynamic Qt5
    set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ -lwsock32 -lws2_32 ${CMAKE_CXX_STANDARD_LIBRARIES}")
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive,-Bdynamic")
    set(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc -static-libstdc++")

    set(USE_STAT OFF CACHE BOOL "Disable stat() for Windows" FORCE)
    # Enable readline for Windows since it's available in MinGW environment
    set(USE_READLINE ON CACHE BOOL "Enable readline for Windows" FORCE)
    add_definitions(-DWINDOWS_BUILD)
    
    # Set crosscompiling flag for Minuit2
    set(CMAKE_CROSSCOMPILING TRUE)
    message(STATUS "Configuring for Windows cross-compilation with static runtime linking")
endif()

# GSL setup - different approach for Windows cross-compilation
if(CROSS_COMPILE_WINDOWS)
    # For Windows cross-compilation, find static libraries explicitly
    set(GSL_INCLUDE_DIRS "/usr/x86_64-w64-mingw32/sys-root/mingw/include")
    
    # Try to find static libraries first, then fall back to dynamic
    find_library(GSL_STATIC_LIBRARY 
        NAMES libgsl.a
        PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib"
        NO_DEFAULT_PATH)
    
    find_library(GSLCBLAS_STATIC_LIBRARY 
        NAMES libgslcblas.a
        PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib"
        NO_DEFAULT_PATH)
    
    if(GSL_STATIC_LIBRARY AND GSLCBLAS_STATIC_LIBRARY)
        set(GSL_LIBRARIES ${GSL_STATIC_LIBRARY} ${GSLCBLAS_STATIC_LIBRARY} m)
        message(STATUS "Found GSL static libraries: ${GSL_LIBRARIES}")
    else()
        # Fall back to dynamic libraries if static not available
        find_library(GSL_DYNAMIC_LIBRARY 
            NAMES libgsl.dll.a gsl
            PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib"
            NO_DEFAULT_PATH)
        
        find_library(GSLCBLAS_DYNAMIC_LIBRARY 
            NAMES libgslcblas.dll.a gslcblas
            PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib"
            NO_DEFAULT_PATH)
        
        if(GSL_DYNAMIC_LIBRARY AND GSLCBLAS_DYNAMIC_LIBRARY)
            set(GSL_LIBRARIES ${GSL_DYNAMIC_LIBRARY} ${GSLCBLAS_DYNAMIC_LIBRARY} m)
            message(STATUS "Found GSL dynamic libraries: ${GSL_LIBRARIES}")
        else()
            # Final fallback to library names
            set(GSL_LIBRARIES "gsl;gslcblas;m")
            message(WARNING "GSL libraries not found, using library names")
        endif()
    endif()
    
    include_directories(${GSL_INCLUDE_DIRS})
else()
    # Native compilation uses pkg-config
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GSL REQUIRED gsl)
    include_directories(${GSL_INCLUDE_DIRS})
    link_directories(${GSL_LIBRARY_DIRS})
endif()

# OpenMP setup
if(CROSS_COMPILE_WINDOWS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    set(OpenMP_FOUND TRUE)
    message(STATUS "Using OpenMP for Windows cross-compilation")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    # Ask the toolchain first.  GCC, conda's clang and MacPorts all answer here,
    # and whatever they provide is the runtime their standard library expects.
    find_package(OpenMP QUIET)
    if(NOT OpenMP_CXX_FOUND AND APPLE AND _libomp_prefix)
        # Only Apple's clang gets this far: it understands no -fopenmp, and its
        # runtime is a separate Homebrew formula.  Supply what FindOpenMP cannot
        # infer, then ask again.
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${_libomp_prefix}/include")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${_libomp_prefix}/include")
        set(OpenMP_C_LIB_NAMES omp)
        set(OpenMP_CXX_LIB_NAMES omp)
        set(OpenMP_omp_LIBRARY "${_libomp_prefix}/lib/libomp.dylib")
        message(STATUS "OpenMP: falling back to Homebrew libomp at ${_libomp_prefix}")
        find_package(OpenMP QUIET)
    endif()
    if(NOT OpenMP_CXX_FOUND)
        message(FATAL_ERROR
            "OpenMP not found.  On macOS install it with 'brew install libomp' "
            "(Apple's clang ships no OpenMP runtime); on Debian/Ubuntu it comes "
            "with build-essential.")
    endif()
    message(STATUS "Using OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

    # Finding OpenMP is not the same as being able to include <omp.h>, and seven
    # source files include it.  For AppleClang, FindOpenMP locates the library
    # and reports "-Xclang -fopenmp" with no -I and an empty
    # OpenMP_CXX_INCLUDE_DIRS, while Homebrew keeps libomp keg-only -- so the
    # header sits in a directory nothing searches and the first compile fails
    # with "'omp.h' file not found" despite a successful configure.
    #
    # Test rather than assume: a conda or MacPorts toolchain has its own omp.h
    # already on the include path and must keep using it, since its header has
    # to match the runtime its standard library links.
    include(CheckIncludeFileCXX)
    get_property(_azure2_incs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
    set(CMAKE_REQUIRED_FLAGS "${OpenMP_CXX_FLAGS}")
    set(CMAKE_REQUIRED_INCLUDES ${_azure2_incs} ${OpenMP_CXX_INCLUDE_DIRS})
    check_include_file_cxx(omp.h AZURE2_HAVE_OMP_H)
    unset(CMAKE_REQUIRED_FLAGS)
    unset(CMAKE_REQUIRED_INCLUDES)

    if(OpenMP_CXX_INCLUDE_DIRS)
        include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    endif()
    if(NOT AZURE2_HAVE_OMP_H)
        if(_libomp_prefix AND EXISTS "${_libomp_prefix}/include/omp.h")
            include_directories("${_libomp_prefix}/include")
            message(STATUS "OpenMP: omp.h taken from ${_libomp_prefix}/include")
        else()
            message(FATAL_ERROR
                "OpenMP was found but <omp.h> is not on the include path.  On "
                "macOS install the runtime with 'brew install libomp'; if it is "
                "already installed somewhere unusual, pass its location with "
                "-DCMAKE_PREFIX_PATH=/path/to/libomp.")
        endif()
    endif()
endif()

# Extra includes/libs
set(AZURE2_SOURCE_DIR ${CMAKE_BUILD_DIRECTORY})
include_directories(${AZURE2_SOURCE_DIR}/external/coul/include)
# OpenMP_LIBRARIES is not a variable FindOpenMP sets -- the per-language
# OpenMP_CXX_LIBRARIES is -- so this list silently carried nothing.  That went
# unnoticed wherever the compiler driver links the runtime itself from
# -fopenmp (GCC, and clang outside Apple's), but AppleClang is driven with
# -Xpreprocessor -fopenmp, which is a preprocessor flag only: nothing links
# libomp, and the build fails at the link step on undefined _omp_* symbols.
set(EXTRA_LIBS coul pugixml ${OpenMP_CXX_LIBRARIES})

# Windows DLL bundling function (define before subdirectories)
function(copy_dll_if_exists target_name dll_name dll_path)
    if(EXISTS "${dll_path}")
        add_custom_command(TARGET ${target_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E echo "=== COPYING DLL: ${dll_name} ==="
            COMMAND ${CMAKE_COMMAND} -E echo "FROM: ${dll_path}"
            COMMAND ${CMAKE_COMMAND} -E echo "TO: $<TARGET_FILE_DIR:${target_name}>/${dll_name}"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${dll_path}" $<TARGET_FILE_DIR:${target_name}>
            COMMAND ${CMAKE_COMMAND} -E echo "=== DLL COPY COMPLETED: ${dll_name} ==="
            COMMENT "Copying ${dll_name}")
    else()
        message(WARNING "DLL not found: ${dll_path}")
    endif()
endfunction()

# Minuit2 build
add_subdirectory(external/minuit2)

# Add Minuit2 to include directories
include_directories(${AZURE2_SOURCE_DIR}/external/minuit2/inc)

# NLopt build (optional minimizer)
option(USE_NLOPT "Use NLopt library for minimization" OFF)
if(USE_NLOPT)
    # Configure NLopt to build as static library with minimal features
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(NLOPT_CXX ON CACHE BOOL "" FORCE)
    set(NLOPT_FORTRAN OFF CACHE BOOL "" FORCE)
    set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
    set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
    set(NLOPT_MATLAB OFF CACHE BOOL "" FORCE)
    set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
    set(NLOPT_JAVA OFF CACHE BOOL "" FORCE)
    set(NLOPT_SWIG OFF CACHE BOOL "" FORCE)
    set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)

    add_subdirectory(nlopt)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/nlopt/src/api)
    add_definitions(-DUSE_NLOPT)
    message(STATUS "NLopt minimization enabled")
else()
    message(STATUS "NLopt minimization disabled - using Minuit2 only")
endif()

# ERYA SRIM utilities (optional)
if(USE_ERYA)
    # Add ERYA subdirectory (pugixml is handled inside external/erya/CMakeLists.txt)
    add_subdirectory(external/erya)
    
    # Add include directories from ERYA module
    include_directories(${SRIM_UTILITIES_INCLUDE_DIRS})
    
    # Define USE_ERYA preprocessor macro
    add_definitions(-DUSE_ERYA)
    
    message(STATUS "ERYA SRIM utilities enabled")
else()
    message(STATUS "ERYA SRIM utilities disabled")
endif()

# MCMC Bayesian sampling (optional)
if(USE_MCMC)
    # Check if numcmc library directory exists
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/numcmc/CMakeLists.txt")
        # Add numcmc as subdirectory
        add_subdirectory(numcmc)
        
        # Set up include directories
        include_directories("${CMAKE_CURRENT_SOURCE_DIR}/numcmc/include")
        
        # Define USE_MCMC preprocessor macro
        add_definitions(-DUSE_MCMC)
        
        # Add numcmc to extra libs for linking
        set(EXTRA_LIBS ${EXTRA_LIBS} numcmc)
        
        # Set RPATH for the executable to find libnumcmc.so at runtime
        set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
        set(CMAKE_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR}/numcmc")
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
        
        message(STATUS "MCMC Bayesian sampling enabled with numcmc library")
    else()
        message(FATAL_ERROR "numcmc library not found. Please clone the numcmc library to the numcmc/ directory")
    endif()
else()
    message(STATUS "MCMC Bayesian sampling disabled")
endif()

# GUI and core modules
if(USE_API)
    add_subdirectory(api)
endif()

# Qt5 setup for cross-compilation
if(BUILD_GUI)
    if(CROSS_COMPILE_WINDOWS)
        # For Windows cross-compilation, set Qt5 paths manually
        set(Qt5_DIR "/usr/x86_64-w64-mingw32/sys-root/mingw/lib/cmake/Qt5")
        set(CMAKE_PREFIX_PATH "/usr/x86_64-w64-mingw32/sys-root/mingw")
        
        # Add explicit paths for Qt5 dependency libraries to help linker
        link_directories("/usr/x86_64-w64-mingw32/sys-root/mingw/lib")
        
        # Find Qt5 components needed for Windows
        find_package(Qt5 REQUIRED COMPONENTS Core Widgets)
        
        # Explicitly link Qt5 dependency libraries to ensure they're found
        # This helps with the DLL loading issue by ensuring import libraries are linked
        find_library(ICONV_LIB NAMES iconv PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib" NO_DEFAULT_PATH)
        find_library(PCRE2_LIB NAMES pcre2-16 PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib" NO_DEFAULT_PATH)
        find_library(ZLIB_LIB NAMES z PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib" NO_DEFAULT_PATH)
        find_library(HARFBUZZ_LIB NAMES harfbuzz PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib" NO_DEFAULT_PATH)
        find_library(PNG_LIB NAMES png16 PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib" NO_DEFAULT_PATH)
        
        if(ICONV_LIB)
            message(STATUS "Found iconv library: ${ICONV_LIB}")
            set(EXTRA_LIBS ${EXTRA_LIBS} ${ICONV_LIB})
        endif()
        if(PCRE2_LIB)
            message(STATUS "Found pcre2-16 library: ${PCRE2_LIB}")
            set(EXTRA_LIBS ${EXTRA_LIBS} ${PCRE2_LIB})
        endif()
        if(ZLIB_LIB)
            message(STATUS "Found zlib library: ${ZLIB_LIB}")
            set(EXTRA_LIBS ${EXTRA_LIBS} ${ZLIB_LIB})
        endif()
        if(HARFBUZZ_LIB)
            message(STATUS "Found harfbuzz library: ${HARFBUZZ_LIB}")
            set(EXTRA_LIBS ${EXTRA_LIBS} ${HARFBUZZ_LIB})
        endif()
        if(PNG_LIB)
            message(STATUS "Found png library: ${PNG_LIB}")
            set(EXTRA_LIBS ${EXTRA_LIBS} ${PNG_LIB})
        endif()
        
        # Handle QWT for Windows if enabled
        if(USE_QWT)
            set(QWT_INCLUDE_DIRS "/usr/x86_64-w64-mingw32/sys-root/mingw/include/qt5/qwt")
            set(QWT_LIBRARIES "qwt-qt5")
        endif()
    else()
        # Native Linux/macOS QWT configuration
        if(USE_QWT)
            # Set the CMake module path to find our custom FindQwt.cmake
            list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
            
            # Find QWT using our custom module
            find_package(Qwt REQUIRED)
            
            if(Qwt_FOUND)
                message(STATUS "Found QWT: ${QWT_LIBRARIES}")
                message(STATUS "QWT include directories: ${QWT_INCLUDE_DIR}")
                
                # Set variables for use in subdirectories
                set(QWT_INCLUDE_DIRS ${QWT_INCLUDE_DIR})
                set(QWT_LIBRARIES ${QWT_LIBRARY})
                
                # Add to global include directories
                include_directories(${QWT_INCLUDE_DIR})
            else()
                message(FATAL_ERROR "QWT library not found. Please install libqwt-qt5-dev or disable USE_QWT")
            endif()
        endif()
    endif()
    
    add_subdirectory(gui)
    add_definitions(-DGUI_BUILD)
endif()

add_subdirectory(external/coul)
add_subdirectory(src)

# Headless GUI tests. After src/, because they link the same libraries the
# AZURE2 executable does; the subdirectory returns immediately unless the GUI
# is being built.
enable_testing()
add_subdirectory(tests/gui)
add_subdirectory(tests/pyazr)

# Install and copy SRIM data file (only if ERYA is enabled)
if(USE_ERYA)
    install(FILES external/erya/data/SRIM2013.xml DESTINATION share/azure2/ COMPONENT Runtime)
    
    # Copy SRIM data file to build directory for development
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/external/erya/data/SRIM2013.xml
                   ${CMAKE_CURRENT_BINARY_DIR}/SRIM2013.xml COPYONLY)
endif()

# Conditional defines
if(NOT USE_STAT)
    add_definitions(-DNO_STAT)
endif()

if(NOT USE_READLINE)
    add_definitions(-DNO_READLINE)
    set(NO_READLINE "TRUE")
else()
    set(NO_READLINE "FALSE")
endif()

# Code coverage
if(USE_GCOV)
    set(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
    set(GCC_COVERAGE_LINK_FLAGS "-lgcov")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
endif()

# Packaging
include(InstallRequiredSystemLibraries)

if(BUILD_MACOS_BUNDLE AND APPLE)
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_DMG_FORMAT "UDBZ")
    set(CPACK_DMG_VOLUME_NAME "AZURE2")
    set(CPACK_SYSTEM_NAME "macOS")
    set(MACOSX_BUNDLE_ICON_FILE azure2.icns)
    set(AZURE2_ICON_PATH "${CMAKE_CURRENT_SOURCE_DIR}/gui/images/azure2.icns")
    install(FILES ${AZURE2_ICON_PATH} DESTINATION "AZURE2.app/Contents/Resources" COMPONENT Runtime)

    set(AZURE2_BUNDLE_PROPERTIES
        MACOSX_BUNDLE TRUE
        MACOSX_BUNDLE_INFO_STRING "AZURE2 R-Matrix Analysis Software"
        MACOSX_BUNDLE_BUNDLE_NAME "AZURE2"
        MACOSX_BUNDLE_BUNDLE_VERSION "${AZURE2_VERSION}"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "${AZURE2_VERSION}"
        MACOSX_BUNDLE_LONG_VERSION_STRING "Version ${AZURE2_VERSION}"
        MACOSX_BUNDLE_ICON_FILE "azure2.icns"
        MACOSX_BUNDLE_GUI_IDENTIFIER "org.azure2.rmatrix"
        MACOSX_BUNDLE_COPYRIGHT "University of Notre Dame")
endif()

if(CROSS_COMPILE_WINDOWS)
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_SYSTEM_NAME "Windows")
    set(CPACK_NSIS_DISPLAY_NAME "AZURE2 R-Matrix Analysis")
    set(CPACK_NSIS_PACKAGE_NAME "AZURE2")
    set(CPACK_NSIS_HELP_LINK "https://azure.nd.edu/")
    set(CPACK_NSIS_URL_INFO_ABOUT "https://azure.nd.edu/")
    set(CPACK_NSIS_CONTACT "azure@nd.edu")
    set(CPACK_NSIS_MODIFY_PATH ON)
    set_target_properties(AZURE2 PROPERTIES SUFFIX ".exe")
    
    # No DLL installation needed for static executable
    # All libraries are embedded in the .exe file
endif()

# General packaging info
set(CPACK_PACKAGE_NAME "AZURE2")
set(CPACK_PACKAGE_VENDOR "University of Notre Dame")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "AZURE2 R-Matrix Analysis Software")
set(CPACK_PACKAGE_VERSION "${AZURE2_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${AZURE2_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${AZURE2_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${AZURE2_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "AZURE2")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

if(NOT CPACK_GENERATOR)
    if(UNIX AND NOT APPLE)
        set(CPACK_GENERATOR "TGZ;DEB")
    endif()
endif()

include(CPack)
