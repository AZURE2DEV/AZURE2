set(AZURE2_SOURCES AChannel.cpp CoulFunc.cpp GenMatrixFunc.cpp
KLGroup.cpp ALevel.cpp Decay.cpp AngCoeff.cpp MGroup.cpp
AMatrixFunc.cpp  DoubleFactorial.cpp EigenFunc.cpp  NFIntegral.cpp
AZURE2.cpp ECIntegral.cpp ECIntegralCache.cpp ECAmplitudeCache.cpp MatrixInv.cpp PPair.cpp AZURECalc.cpp
ECMGroup.cpp ReactionRate.cpp AZUREMain.cpp EData.cpp ShftFunc.cpp
RMatrixFunc.cpp AZUREOutput.cpp  EDataIterator.cpp   TargetEffect.cpp
AZUREParams.cpp EPoint.cpp Interference.cpp CNuc.cpp Equation.cpp
JGroup.cpp Config.cpp ESegment.cpp KGroup.cpp IntegratedFermiFunc.cpp
EffectiveCharge.cpp  GSLException.cpp ParameterLimitsManager.cpp)

include_directories(${AZURE2_SOURCE_DIR}/include)

#if (APPLE)
#set(AZURE2_PROPERTIES MACOSX_BUNDLE)
#endif (APPLE)

if(NOT NO_READLINE)
set (EXTRA_LIBS ${EXTRA_LIBS} readline)
endif (NOT NO_READLINE)
find_library(NCURSES_LIBRARY NAMES ncurses)
if (NCURSES_LIBRARY)
set (EXTRA_LIBS ${EXTRA_LIBS} ${NCURSES_LIBRARY})
mark_as_advanced (NCURSES_LIBRARY)
endif (NCURSES_LIBRARY)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
if (NOT NO_READLINE)
set(EXTRA_LIBS ${EXTRA_LIBS} termcap lcms)
endif (NOT NO_READLINE)
set(AZURE2_SOURCES ${AZURE2_SOURCES} AZURE2.rc)
set(CMAKE_EXE_LINKER_FLAGS "-Wl,-subsystem,windows")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if (BUILD_LIBRARY) 
add_library(AZURE2 ${AZURE2_SOURCES})
else (BUILD_LIBRARY)
add_executable (AZURE2 ${AZURE2_PROPERTIES} ${AZURE2_SOURCES})
endif (BUILD_LIBRARY)

target_link_libraries (AZURE2 ${GSL_LIBRARIES} ${MINUIT2_LIBRARIES} ${AZURESETUP_LIBRARIES} ${AZUREAPI_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBS})

install(TARGETS AZURE2 DESTINATION ${CMAKE_INSTALL_DIR})
