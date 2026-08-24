# Core sources, shared by the AZURE2 executable and the _azure2 Python module.
set(AZURE2_CORE_SOURCES AChannel.cpp CoulFunc.cpp CoulFuncCache.cpp GenMatrixFunc.cpp
KLGroup.cpp ALevel.cpp Decay.cpp AngCoeff.cpp MGroup.cpp
AMatrixFunc.cpp  DoubleFactorial.cpp EigenFunc.cpp  NFIntegral.cpp
ECIntegral.cpp ECAmplitudeCache.cpp MatrixInv.cpp LevelMatrixSolver.cpp
PPair.cpp AZURECalc.cpp
ECMGroup.cpp ReactionRate.cpp AZUREMain.cpp EData.cpp ShftFunc.cpp
RMatrixFunc.cpp AZUREOutput.cpp  EDataIterator.cpp   TargetEffect.cpp
AZUREParams.cpp EPoint.cpp Interference.cpp CNuc.cpp Equation.cpp
JGroup.cpp Config.cpp ESegment.cpp KGroup.cpp IntegratedFermiFunc.cpp
EffectiveCharge.cpp  GSLException.cpp ParameterLimitsManager.cpp Straggling.cpp
AdaptiveIntegrationGrid.cpp NuclearPotentialManager.cpp AZUREGrad.cpp ThmFunc.cpp
THMMatrixFunc.cpp CovarianceBand.cpp ParameterLabel.cpp PolarizationFunc.cpp)

# Add GUI stub when GUI is not enabled
if(NOT BUILD_GUI)
    list(APPEND AZURE2_CORE_SOURCES gui_stub.cpp)
endif()

# Add MCMC sources when MCMC is enabled
if(USE_MCMC)
    list(APPEND AZURE2_CORE_SOURCES AZURECalcMCMC.cpp)
endif()

include_directories(${AZURE2_SOURCE_DIR}/include)

# The engine as a library, so both the CLI/GUI binary and the _azure2 Python
# module link the same compiled code.  POSITION_INDEPENDENT_CODE lets the
# static objects be pulled into the Python extension.
add_library(AZURE2Core STATIC ${AZURE2_CORE_SOURCES})
target_include_directories(AZURE2Core PUBLIC ${AZURE2_SOURCE_DIR}/include)
set_target_properties(AZURE2Core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# macOS Bundle properties (set from parent CMakeLists.txt)
if(BUILD_MACOS_BUNDLE AND APPLE)
    set(AZURE2_PROPERTIES MACOSX_BUNDLE)
endif()

if(NOT NO_READLINE)
    if(CROSS_COMPILE_WINDOWS)
        # For Windows cross-compilation, link directly to the MinGW readline library files
        find_library(READLINE_LIBRARY 
            NAMES libreadline.dll.a
            PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib"
            NO_DEFAULT_PATH)
        
        find_library(HISTORY_LIBRARY 
            NAMES libhistory.dll.a
            PATHS "/usr/x86_64-w64-mingw32/sys-root/mingw/lib"
            NO_DEFAULT_PATH)
        
        if(READLINE_LIBRARY)
            set (EXTRA_LIBS ${EXTRA_LIBS} ${READLINE_LIBRARY})
            message(STATUS "Found readline library: ${READLINE_LIBRARY}")
        endif()
        
        if(HISTORY_LIBRARY)
            set (EXTRA_LIBS ${EXTRA_LIBS} ${HISTORY_LIBRARY})
            message(STATUS "Found history library: ${HISTORY_LIBRARY}")
        endif()
    else()
        # For native compilation, use standard readline
        set (EXTRA_LIBS ${EXTRA_LIBS} readline)
        
        find_library(NCURSES_LIBRARY NAMES ncurses)
        if (NCURSES_LIBRARY)
            set (EXTRA_LIBS ${EXTRA_LIBS} ${NCURSES_LIBRARY})
            mark_as_advanced (NCURSES_LIBRARY)
        endif (NCURSES_LIBRARY)
    endif()
endif (NOT NO_READLINE)

set(AZURE2_EXEC_SOURCES AZURE2.cpp)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(AZURE2_EXEC_SOURCES ${AZURE2_EXEC_SOURCES} AZURE2.rc)
    # Note: For console applications, use -subsystem,console or remove subsystem flag entirely
    # We'll let the static linking flags from main CMakeLists.txt handle the linker configuration
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if (BUILD_LIBRARY) 
add_library(AZURE2 STATIC ${AZURE2_CORE_SOURCES} ${AZURE2_EXEC_SOURCES})
else (BUILD_LIBRARY)
add_executable (AZURE2 ${AZURE2_PROPERTIES} ${AZURE2_EXEC_SOURCES})
endif (BUILD_LIBRARY)

# The engine's own dependencies, declared on the engine and PUBLIC so that
# everything linking it -- the executable, the GUI, the Python module -- gets
# them placed *after* it automatically.
#
# This has to be derived rather than hand-written.  A static archive only
# supplies the symbols already demanded when the linker reaches it, so listing
# Minuit2 or GSL before AZURE2Core leaves the engine's references to them
# unresolved on GNU ld (Linux and MinGW).  macOS's linker rescans archives and
# hides the mistake, which is why the hand-written order looked fine.
target_link_libraries(AZURE2Core PUBLIC
    Minuit2::Minuit2
    ${GSL_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${EXTRA_LIBS})
if(USE_NLOPT)
    target_link_libraries(AZURE2Core PUBLIC nlopt)
endif()

# AZURESetup (the GUI) calls into the engine, so it precedes it; AZUREAPI in
# turn pulls AZURE2Core in as a PUBLIC dependency of its own.  Everything else
# follows from the graph above.
target_link_libraries(AZURE2 ${AZURESETUP_LIBRARIES} ${AZUREAPI_LIBRARIES} AZURE2Core)

# Apply macOS bundle properties if building bundle
if(BUILD_MACOS_BUNDLE AND APPLE AND DEFINED AZURE2_BUNDLE_PROPERTIES)
    set_target_properties(AZURE2 PROPERTIES ${AZURE2_BUNDLE_PROPERTIES})
    
    # Set RPATH for better dependency management
    set_target_properties(AZURE2 PROPERTIES
        MACOSX_RPATH TRUE
        BUILD_WITH_INSTALL_RPATH FALSE
        INSTALL_RPATH "@loader_path/../Libraries;@loader_path/../Frameworks;@executable_path/../Libraries;@executable_path/../Frameworks"
        INSTALL_RPATH_USE_LINK_PATH TRUE)
    
    # Copy icon file to Resources during build
    if(EXISTS "${AZURE2_ICON_PATH}")
        add_custom_command(TARGET AZURE2 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_BUNDLE_DIR:AZURE2>/Contents/Resources"
            COMMAND ${CMAKE_COMMAND} -E copy "${AZURE2_ICON_PATH}" "$<TARGET_BUNDLE_DIR:AZURE2>/Contents/Resources/"
            COMMENT "Copying icon to bundle Resources")
    endif()
    
    # Create Libraries directory for bundled dependencies
    add_custom_command(TARGET AZURE2 POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_BUNDLE_DIR:AZURE2>/Contents/Libraries"
        COMMENT "Creating Libraries directory for dependencies")
endif()

# Windows DLL bundling for libraries that can't be statically linked
if(CROSS_COMPILE_WINDOWS)
    set(MINGW_DLL_PATH "/usr/x86_64-w64-mingw32/sys-root/mingw/bin")
    
    # Add application manifest for better DLL loading (MinGW approach)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/AZURE2.manifest")
        # Copy manifest to output directory
        add_custom_command(TARGET AZURE2 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/AZURE2.manifest"
            "$<TARGET_FILE_DIR:AZURE2>/AZURE2.exe.manifest"
            COMMENT "Copying application manifest")
    endif()
    
    # Copy essential MinGW runtime DLLs (fallback in case static linking doesn't work)
    copy_dll_if_exists(AZURE2 "libgcc_s_seh-1.dll" "${MINGW_DLL_PATH}/libgcc_s_seh-1.dll")
    copy_dll_if_exists(AZURE2 "libstdc++-6.dll" "${MINGW_DLL_PATH}/libstdc++-6.dll")
    copy_dll_if_exists(AZURE2 "libwinpthread-1.dll" "${MINGW_DLL_PATH}/libwinpthread-1.dll")
    
    # Copy DLLs for libraries that must be dynamically linked
    # GSL DLL (if static libraries not found)
    copy_dll_if_exists(AZURE2 "libgsl-0.dll" "${MINGW_DLL_PATH}/libgsl-0.dll")
    copy_dll_if_exists(AZURE2 "libgslcblas-0.dll" "${MINGW_DLL_PATH}/libgslcblas-0.dll")
    
    # Readline DLLs (these are dynamically linked via .dll.a)
    copy_dll_if_exists(AZURE2 "libreadline8.dll" "${MINGW_DLL_PATH}/libreadline8.dll")
    copy_dll_if_exists(AZURE2 "libhistory8.dll" "${MINGW_DLL_PATH}/libhistory8.dll")
    
    # Qt5 DLLs (these are always dynamically linked in MinGW)
    if(BUILD_GUI)
        copy_dll_if_exists(AZURE2 "Qt5Core.dll" "${MINGW_DLL_PATH}/Qt5Core.dll")
        copy_dll_if_exists(AZURE2 "Qt5Gui.dll" "${MINGW_DLL_PATH}/Qt5Gui.dll")
        copy_dll_if_exists(AZURE2 "Qt5Widgets.dll" "${MINGW_DLL_PATH}/Qt5Widgets.dll")
        copy_dll_if_exists(AZURE2 "Qt5PrintSupport.dll" "${MINGW_DLL_PATH}/Qt5PrintSupport.dll")
        copy_dll_if_exists(AZURE2 "Qt5OpenGL.dll" "${MINGW_DLL_PATH}/Qt5OpenGL.dll")
        copy_dll_if_exists(AZURE2 "Qt5Svg.dll" "${MINGW_DLL_PATH}/Qt5Svg.dll")
        
        # Qt5 dependency DLLs (required by Qt5Core.dll and Qt5Gui.dll)
        copy_dll_if_exists(AZURE2 "iconv.dll" "${MINGW_DLL_PATH}/iconv.dll")
        copy_dll_if_exists(AZURE2 "libiconv-2.dll" "${MINGW_DLL_PATH}/libiconv-2.dll")
        copy_dll_if_exists(AZURE2 "libpcre2-16-0.dll" "${MINGW_DLL_PATH}/libpcre2-16-0.dll")
        copy_dll_if_exists(AZURE2 "libpcre2-8-0.dll" "${MINGW_DLL_PATH}/libpcre2-8-0.dll")
        copy_dll_if_exists(AZURE2 "zlib1.dll" "${MINGW_DLL_PATH}/zlib1.dll")
        copy_dll_if_exists(AZURE2 "libharfbuzz-0.dll" "${MINGW_DLL_PATH}/libharfbuzz-0.dll")
        copy_dll_if_exists(AZURE2 "libpng16-16.dll" "${MINGW_DLL_PATH}/libpng16-16.dll")
        copy_dll_if_exists(AZURE2 "libfreetype-6.dll" "${MINGW_DLL_PATH}/libfreetype-6.dll")
        copy_dll_if_exists(AZURE2 "libglib-2.0-0.dll" "${MINGW_DLL_PATH}/libglib-2.0-0.dll")
        copy_dll_if_exists(AZURE2 "libintl-8.dll" "${MINGW_DLL_PATH}/libintl-8.dll")
        copy_dll_if_exists(AZURE2 "libgraphite2.dll" "${MINGW_DLL_PATH}/libgraphite2.dll")
        copy_dll_if_exists(AZURE2 "libbz2-1.dll" "${MINGW_DLL_PATH}/libbz2-1.dll")
        
        # Additional readline dependency DLL
        copy_dll_if_exists(AZURE2 "libtermcap-0.dll" "${MINGW_DLL_PATH}/libtermcap-0.dll")
        
        # QWT DLL
        if(USE_QWT)
            copy_dll_if_exists(AZURE2 "qwt-qt5.dll" "${MINGW_DLL_PATH}/qwt-qt5.dll")
        endif()
        
        # Qt5 platform plugins (required for Qt GUI)
        add_custom_command(TARGET AZURE2 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E echo "=== COPYING QT5 PLATFORM PLUGIN ==="
            COMMAND ${CMAKE_COMMAND} -E echo "FROM: ${MINGW_DLL_PATH}/../lib/qt5/plugins/platforms/qwindows.dll"
            COMMAND ${CMAKE_COMMAND} -E echo "TO: $<TARGET_FILE_DIR:AZURE2>/platforms/"
            COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:AZURE2>/platforms"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${MINGW_DLL_PATH}/../lib/qt5/plugins/platforms/qwindows.dll"
            "$<TARGET_FILE_DIR:AZURE2>/platforms/"
            COMMAND ${CMAKE_COMMAND} -E echo "=== QT5 PLATFORM PLUGIN COPY COMPLETED ==="
            COMMENT "Copying Qt5 platform plugin")
    endif()
    
    add_custom_command(TARGET AZURE2 POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E echo "=== WINDOWS BUILD SUMMARY ==="
        COMMAND ${CMAKE_COMMAND} -E echo "TARGET_FILE_DIR for AZURE2: $<TARGET_FILE_DIR:AZURE2>"
        COMMAND ${CMAKE_COMMAND} -E echo "Contents of target directory:"
        COMMAND ${CMAKE_COMMAND} -E echo "Listing files in $<TARGET_FILE_DIR:AZURE2>:"
        COMMAND ls -la "$<TARGET_FILE_DIR:AZURE2>/" || echo "Could not list target directory"
        COMMAND ${CMAKE_COMMAND} -E echo "=== END BUILD SUMMARY ==="
        COMMENT "Windows executable with DLLs built")
endif()

# Install target with proper bundle handling
if(BUILD_MACOS_BUNDLE AND APPLE)
    install(TARGETS AZURE2 
            BUNDLE DESTINATION "."
            RUNTIME DESTINATION "bin")
else()
    install(TARGETS AZURE2 DESTINATION ${CMAKE_INSTALL_DIR})
endif()
