# physics/: the physical objects and stages of the R-matrix calculation (docs/upgrade/ARCHITECTURE.md sections 2-3).
# It uses nothing else of AZURE2: no experiment/, inference/, io/, engine (src/), Qt, GSL or OpenMP.  The rule is
# enforced three times: the target links nothing and sees only its own include directory; the ctest
# physics_standalone configures and builds this directory as a project of its own (tests/physics_standalone.cmake),
# where no AZURE2 header or library exists to be found; and repo_consistency (tests/consistency) refuses any
# include but azure2/physics/ headers and the C++ standard library, which also catches a GSL, Qt or OpenMP header
# that a standalone build would find on the system include path.
#
# Layout: public headers in include/azure2/physics/ (included as "azure2/physics/X.h"), sources in src/, tests that
# need only physics/ in tests/.
#
# It works as a subdirectory of AZURE2 and as a top-level project (cmake -S physics -B build-physics).
cmake_minimum_required(VERSION 3.16...4.0)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(azure2_physics LANGUAGES CXX)
  set(AZURE2_PHYSICS_STANDALONE ON)
  enable_testing()
endif()

add_library(azure2_physics STATIC
  src/AngularCoupling.cpp
  src/JpiBlock.cpp
  src/SymmetricEigen.cpp
  src/System.cpp)
add_library(azure2::physics ALIAS azure2_physics)
target_include_directories(azure2_physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(azure2_physics PUBLIC cxx_std_17)
set_target_properties(azure2_physics PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_physics PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Unit test of the types themselves (construction, validation, index conventions); links physics only.
add_executable(physics_system_unit_test tests/system_unit_test.cpp)
target_link_libraries(physics_system_unit_test PRIVATE azure2_physics)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(physics_system_unit_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME physics_system_unit COMMAND physics_system_unit_test)

# The eigensolvers against their defining identities.
# (Sources named in full: repo_consistency checks that every tests/*_test.cpp is named here.)
foreach(_src symmetric_eigen_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(physics_${_t}_test tests/${_src})
  target_link_libraries(physics_${_t}_test PRIVATE azure2_physics)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(physics_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME physics_${_t} COMMAND physics_${_t}_test)
endforeach()

if(NOT AZURE2_PHYSICS_STANDALONE)
  # physics/ configured, built and tested on its own, in a fresh build tree under this one.
  add_test(NAME physics_standalone
           COMMAND ${CMAKE_COMMAND}
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}
                   -DBIN=${CMAKE_CURRENT_BINARY_DIR}/standalone
                   -DGENERATOR=${CMAKE_GENERATOR}
                   -DCXX=${CMAKE_CXX_COMPILER}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/physics_standalone.cmake)
  set_tests_properties(physics_standalone PROPERTIES TIMEOUT 300)
endif()
