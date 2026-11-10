# Core sources, shared by the AZURE2 executable and the _azure2 Python module.
set(AZURE2_CORE_SOURCES AChannel.cpp CoulFunc.cpp CoulombWavesAt.cpp ChannelTable.cpp InputFiles.cpp GenMatrixFunc.cpp
KLGroup.cpp ALevel.cpp Decay.cpp AngCoeff.cpp MGroup.cpp
AMatrixFunc.cpp  DoubleFactorial.cpp EigenFunc.cpp  NFIntegral.cpp
ECIntegral.cpp ECTable.cpp ECAmplitudeCache.cpp MatrixInv.cpp LevelMatrixSolver.cpp Profile.cpp ChannelSpace.cpp
PPair.cpp AZURECalc.cpp Objective.cpp
ECMGroup.cpp ReactionRate.cpp AZUREMain.cpp EData.cpp DataCovariance.cpp ShftFunc.cpp ScatteringLength.cpp ModelData.cpp CNucObserved.cpp
RMatrixFunc.cpp AZUREOutput.cpp  EDataIterator.cpp   TargetEffect.cpp
AZUREParams.cpp EPoint.cpp Interference.cpp CNuc.cpp Equation.cpp
JGroup.cpp Config.cpp ESegment.cpp ModelPass.cpp CoulombReference.cpp Predictions.cpp CalculationPlan.cpp KGroup.cpp IntegratedFermiFunc.cpp
EffectiveCharge.cpp  GSLException.cpp ParameterLimitsManager.cpp
AdaptiveIntegrationGrid.cpp NuclearPotentialManager.cpp AZUREGrad.cpp ThmFunc.cpp
EDataOutputFiles.cpp
EDataCheckFiles.cpp
EDataCovariance.cpp
EDataTargetEffects.cpp
THMMatrixFunc.cpp ThmVertexBoundary.cpp ThmNumerics.cpp EDataThm.cpp ConfigThm.cpp ThmExperiment.cpp ThmLineshape.cpp ThmDistortion.cpp ThmOptical.cpp ThmDwVertex.cpp ThmAngular.cpp CovarianceBand.cpp PriorCovariance.cpp ParameterLabel.cpp PolarizationFunc.cpp ChannelFunc.cpp CoulombTable.cpp
EDataThmReports.cpp
CNucSystem.cpp CNucLevelParameters.cpp CNucExternalCapture.cpp CNucReport.cpp CNucChecks.cpp
CNucTransformMessages.cpp ProjectSetup.cpp)

# Add MCMC sources when MCMC is enabled
if(USE_MCMC)
    list(APPEND AZURE2_CORE_SOURCES AZURECalcMCMC.cpp)
endif()

include_directories(${AZURE2_SOURCE_DIR}/include)

# The engine as a library, so both the CLI binary and the _azure2 Python
# module link the same compiled code.  POSITION_INDEPENDENT_CODE lets the
# static objects be pulled into the Python extension.
add_library(AZURE2Core STATIC ${AZURE2_CORE_SOURCES})
target_include_directories(AZURE2Core PUBLIC ${AZURE2_SOURCE_DIR}/include)
set_target_properties(AZURE2Core PROPERTIES POSITION_INDEPENDENT_CODE ON)

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

# ProjectCommands.cpp: --convert, --validate, --schema (they use convert/, linked below)
set(AZURE2_EXEC_SOURCES AZURE2.cpp ProjectCommands.cpp)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(AZURE2_EXEC_SOURCES ${AZURE2_EXEC_SOURCES} AZURE2.rc)
    # Note: For console applications, use -subsystem,console or remove subsystem flag entirely
    # We'll let the static linking flags from main CMakeLists.txt handle the linker configuration
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if (BUILD_LIBRARY) 
add_library(AZURE2 STATIC ${AZURE2_CORE_SOURCES} ${AZURE2_EXEC_SOURCES})
else (BUILD_LIBRARY)
add_executable (AZURE2 ${AZURE2_EXEC_SOURCES})
endif (BUILD_LIBRARY)

# The engine's own dependencies, declared on the engine and PUBLIC so that
# everything linking it -- the executable and the Python module -- gets
# them placed *after* it automatically.
#
# This has to be derived rather than hand-written.  A static archive only
# supplies the symbols already demanded when the linker reaches it, so listing
# Minuit2 or GSL before AZURE2Core leaves the engine's references to them
# unresolved on GNU ld (Linux and MinGW).  macOS's linker rescans archives and
# hides the mistake, which is why the hand-written order looked fine.
# azure2_physics (physics/) is below the engine: the engine may use it, never the reverse.  Only the
# CNuc -> System adapter (CNucSystem.cpp) uses it so far.
# azure2_experiment (experiment/: model expressions, lab -> c.m. frames, the format-3 table of AZURE2 --schema
# and the binding's schema()) uses physics/ and io/'s Record.h only and is below the engine likewise.
# azure2_experiment_reader (the reading path of a format-3 project, ProjectReader.h, and the <levels> tables) links io/.
# azure2_convert_specs (convert/: the old sections as set-up specs, with the old line parsers it is built on) until
# step 3d deletes the old sections; the engine reads an old project only through it (src/ProjectSetup.cpp).
target_link_libraries(AZURE2Core PUBLIC
    azure2_physics
    azure2_io
    azure2_experiment
    azure2_experiment_reader
    azure2_convert_specs
    Minuit2::Minuit2
    ${GSL_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${EXTRA_LIBS})

# Builds, validates and describes the universal Coulomb tables (CoulombTable);
# see docs/source/developer/coulomb_table.rst.
add_executable(azure2-coulomb-table coulomb_table_main.cpp)
target_link_libraries(azure2-coulomb-table PRIVATE AZURE2Core)
if(USE_NLOPT)
    target_link_libraries(AZURE2Core PUBLIC nlopt)
endif()

# AZUREAPI pulls AZURE2Core in as a PUBLIC dependency of its own.  Everything
# else follows from the graph above.
target_link_libraries(AZURE2 ${AZUREAPI_LIBRARIES} AZURE2Core)
# convert/: AZURE2 --convert ... --to 3 and --validate (src/AZURE2.cpp)
target_link_libraries(AZURE2 azure2_convert)

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
    
    # Readline's own dependency
    copy_dll_if_exists(AZURE2 "libtermcap-0.dll" "${MINGW_DLL_PATH}/libtermcap-0.dll")
    
    add_custom_command(TARGET AZURE2 POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E echo "=== WINDOWS BUILD SUMMARY ==="
        COMMAND ${CMAKE_COMMAND} -E echo "TARGET_FILE_DIR for AZURE2: $<TARGET_FILE_DIR:AZURE2>"
        COMMAND ${CMAKE_COMMAND} -E echo "Contents of target directory:"
        COMMAND ${CMAKE_COMMAND} -E echo "Listing files in $<TARGET_FILE_DIR:AZURE2>:"
        COMMAND ls -la "$<TARGET_FILE_DIR:AZURE2>/" || echo "Could not list target directory"
        COMMAND ${CMAKE_COMMAND} -E echo "=== END BUILD SUMMARY ==="
        COMMENT "Windows executable with DLLs built")
endif()

install(TARGETS AZURE2 DESTINATION ${CMAKE_INSTALL_DIR})
