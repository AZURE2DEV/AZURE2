# Core sources, shared by the AZURE2 executable and the _azure2 Python module.
set(AZURE2_CORE_SOURCES AChannel.cpp CoulFunc.cpp CoulombWavesAt.cpp ChannelTable.cpp InputFiles.cpp GenMatrixFunc.cpp
ALevel.cpp LevelInputs.cpp LevelState.cpp LevelParameterLayout.cpp LevelStateFill.cpp ObservedLevelTable.cpp PathwayTable.cpp CompoundSystem.cpp CaptureAyRequests.cpp AngCoeff.cpp
AMatrixFunc.cpp  DoubleFactorial.cpp EigenFunc.cpp  NFIntegral.cpp
ECIntegral.cpp ECTable.cpp ECAmplitudeCache.cpp MatrixInv.cpp LevelMatrixSolver.cpp Profile.cpp ChannelSpace.cpp
PPair.cpp AZURECalc.cpp Objective.cpp
ReactionRate.cpp AZUREMain.cpp RunSetUp.cpp ProgramFit.cpp RunOutputs.cpp
Fitters.cpp FitEvaluations.cpp FitMigrad.cpp FitLeastSquares.cpp FitGslTrustRegion.cpp FitNLopt.cpp FitErrors.cpp
WorkingCopy.cpp DataCovariance.cpp ShftFunc.cpp ScatteringLength.cpp ModelData.cpp LevelTransforms.cpp LevelParameterVector.cpp GammaDecayMap.cpp
RMatrixFunc.cpp AZUREOutput.cpp    TargetEffect.cpp
AZUREParams.cpp PosteriorSupport.cpp Unit.cpp SiteObservables.cpp CompoundSystemBuild.cpp AbundanceSamples.cpp Equation.cpp
JGroup.cpp Config.cpp ProjectOptions.cpp Segment.cpp ModelPass.cpp CoulombReference.cpp ReferenceTableLeaf.cpp Predictions.cpp CalculationPlan.cpp EnergyPass.cpp IntegratedFermiFunc.cpp
EffectiveCharge.cpp  GSLException.cpp ParameterLimitsManager.cpp
AdaptiveIntegrationGrid.cpp AdaptiveExtrapolation.cpp AZUREGrad.cpp GradEnergyTasks.cpp ThmFunc.cpp
DataOutputFiles.cpp CommandLineValues.cpp
DataCheckFiles.cpp
DataCovarianceGroups.cpp
TargetSetup.cpp
TargetStopping.cpp
StoppingScale.cpp
ApparatusParameters.cpp
MeshPlanner.cpp
MeshSlots.cpp
ExternalCaptureSetup.cpp
DataParameters.cpp
HandOver.cpp
TaskScratch.cpp
RunScratch.cpp
RunDiagnostics.cpp
WorkingCopyUnits.cpp NodeTable.cpp DataSet.cpp
PointPass.cpp
PointFold.cpp
TargetLayer.cpp
ChannelTables.cpp
UnitLegendre.cpp
UnitExternalCapture.cpp
SegmentParameters.cpp
ShiftFrames.cpp
ShiftedFold.cpp
THMMatrixFunc.cpp ThmCotangents.cpp ThmVertexBoundary.cpp ThmNumerics.cpp ThmNorm.cpp ConfigThm.cpp ThmBlockEdits.cpp ThmExperiment.cpp ThmKeys.cpp ThmLineshape.cpp ThmDistortion.cpp ThmDwVertex.cpp ThmAngular.cpp CovarianceBand.cpp PriorCovariance.cpp ParameterLabel.cpp PolarizationFunc.cpp ChannelFunc.cpp CoulombTable.cpp
ThmReports.cpp
SystemOfCompound.cpp EngineChannelFunctions.cpp ExternalCaptureAmplitudes.cpp ExternalCaptureFinalStates.cpp LevelReport.cpp CompoundSystemChecks.cpp
LevelInputMessages.cpp ProjectSetup.cpp)

# Add MCMC sources when MCMC is enabled
if(USE_MCMC)
    list(APPEND AZURE2_CORE_SOURCES AZURECalcMCMC.cpp SamplerSupport.cpp)
endif()

include_directories(${AZURE2_SOURCE_DIR}/include)

# The engine as a library, so both the CLI binary and the _azure2 Python
# module link the same compiled code.  POSITION_INDEPENDENT_CODE lets the
# static objects be pulled into the Python extension.
add_library(AZURE2Core STATIC ${AZURE2_CORE_SOURCES})
target_include_directories(AZURE2Core PUBLIC ${AZURE2_SOURCE_DIR}/include)
set_target_properties(AZURE2Core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ProjectCommands.cpp: --convert, --convert-params, --validate, --describe, --schema, --print-predictions (the first
# two use convert/, linked below); RunRequest.cpp: the commands and the options of a run, the command line's reader;
# RunPreparation.cpp: a run request into the Config
set(AZURE2_EXEC_SOURCES AZURE2.cpp ProjectCommands.cpp RunRequest.cpp RunPreparation.cpp)

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
# CompoundSystem -> System adapter (SystemOfCompound.cpp) uses it so far.  azure2_physics_transfer (physics/transfer: distorted waves,
# bound states, the zero-range THM weight, the global optical sets) is what the THM distortions call.
# azure2_experiment (experiment/: model expressions, lab -> c.m. frames, the format-3 table of AZURE2 --schema
# and the binding's schema()) uses physics/ and io/'s Record.h only and is below the engine likewise.
# azure2_experiment_reader (the reading path of a format-3 project, ProjectReader.h, and the <levels> tables) links io/.
target_link_libraries(AZURE2Core PUBLIC
    azure2_support
    azure2_physics
    azure2_physics_transfer
    azure2_stopping
    azure2_io
    azure2_experiment
    azure2_experiment_reader
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
# convert/: the converter of old projects, AZURE2 --convert and --convert-params (src/ProjectCommands.cpp); the engine
# library (AZURE2Core) reads format 3 only and links none of convert/
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
    
    
    add_custom_command(TARGET AZURE2 POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E echo "=== WINDOWS BUILD SUMMARY ==="
        COMMAND ${CMAKE_COMMAND} -E echo "TARGET_FILE_DIR for AZURE2: $<TARGET_FILE_DIR:AZURE2>"
        COMMAND ${CMAKE_COMMAND} -E echo "Contents of target directory:"
        COMMAND ${CMAKE_COMMAND} -E echo "Listing files in $<TARGET_FILE_DIR:AZURE2>:"
        COMMAND ls -la "$<TARGET_FILE_DIR:AZURE2>/" || echo "Could not list target directory"
        COMMAND ${CMAKE_COMMAND} -E echo "=== END BUILD SUMMARY ==="
        COMMENT "Windows executable with DLLs built")
endif()

# A wheel build (scikit-build-core) puts the program into the wheel's scripts, so that pip installs it on the PATH.
if(SKBUILD)
    install(TARGETS AZURE2 RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR} COMPONENT python)
else()
    install(TARGETS AZURE2 DESTINATION ${CMAKE_INSTALL_DIR})
endif()
