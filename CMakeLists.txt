cmake_minimum_required (VERSION 4.0)
project (AZURE2)

set(AZURE2_SOURCE_DIR ${CMAKE_BUILD_DIRECTORY})

if (MINUIT_PATH) 
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${MINUIT_PATH})
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${MINUIT_PATH})
endif (MINUIT_PATH) 

if (GSL_PATH)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${GSL_PATH})
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GSL_PATH})
endif (GSL_PATH) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package (gsl REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
find_package (Minuit2 REQUIRED)
include_directories(${MINUIT2_INCLUDE_DIR})
find_package (OpenMP REQUIRED)
include_directories(${OPENMP_INCLUDE_DIR})
#find_package(BLAS REQUIRED)

option (BUILD_LIBRARY
       "Build AZURE2 as a library" OFF)

option (BUILD_GUI 
       "Build AZURE2 with graphical setup program" YES)
option (USE_STAT
          "Use stat() function to check for directories.  Not
          recommended for windows builds" YES)
option (USE_QWT
          "Use Qwt to build Plot tab." YES)
option (USE_READLINE
		  "Use Readline library for CLI input." YES)
option(USE_GCOV 
       "Create a GCov-enabled build." OFF)
option(CODE_COVERAGE 
       "Enable coverage reporting" ON)

add_subdirectory (api)

if (BUILD_GUI)
add_subdirectory (gui)
add_definitions(-DGUI_BUILD)
endif (BUILD_GUI)

include_directories(${AZURE2_SOURCE_DIR}/coul/include)
add_subdirectory (coul)
set (EXTRA_LIBS coul)

if (NOT USE_STAT) 
add_definitions (-DNO_STAT) 
endif (NOT USE_STAT)

if (NOT USE_READLINE)
set (NO_READLINE "TRUE")
add_definitions (-DNO_READLINE)
else (NOT USE_READLINE)
set (NO_READLINE "FALSE")
endif (NOT USE_READLINE)

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O2        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
    -std=c++11
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

if (USE_GCOV)
  set(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
  set(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
endif()
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

add_subdirectory (src)
