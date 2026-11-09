# experiment/: what was measured and what the model must predict (docs/upgrade/data_model.md section 2.5): the model
# expressions and, from step S6, the lab -> c.m. frames.  It uses physics/ (masses, constants, kinematics) and
# nothing else of AZURE2: no io/, engine (src/), GSL, Qt or OpenMP.  The rule is enforced three times, as for
# physics/: the target links only azure2_physics and sees only its own and physics/'s include directories; the ctest
# experiment_standalone configures and builds this directory (with physics/) as a project of its own; and
# repo_consistency refuses any include but azure2/experiment/, azure2/physics/ and the C++ standard library.
#
# Layout: public headers in include/azure2/experiment/ (included as "azure2/experiment/X.h"), sources in src/, tests
# that need only experiment/ and physics/ in tests/.
cmake_minimum_required(VERSION 3.16...4.0)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(azure2_experiment LANGUAGES C CXX)
  set(AZURE2_EXPERIMENT_STANDALONE ON)
  enable_testing()
  # physics/ is the one library experiment/ uses; standalone it is built here from its own directory.
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../physics ${CMAKE_CURRENT_BINARY_DIR}/physics EXCLUDE_FROM_ALL)
endif()

add_library(azure2_experiment STATIC
  src/Expression.cpp)
add_library(azure2::experiment ALIAS azure2_experiment)
target_include_directories(azure2_experiment PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(azure2_experiment PUBLIC azure2_physics)
target_compile_features(azure2_experiment PUBLIC cxx_std_17)
set_target_properties(azure2_experiment PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_experiment PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Each node against the engine's expression written out in the test, bit for bit, and backward against finite
# differences (docs/upgrade/data_model.md step S8).
# (Sources named in full: repo_consistency checks that every tests/*_test.cpp is named here.)
foreach(_src expression_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(experiment_${_t}_test tests/${_src})
  target_link_libraries(experiment_${_t}_test PRIVATE azure2_experiment)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(experiment_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME experiment_${_t} COMMAND experiment_${_t}_test)
endforeach()

if(NOT AZURE2_EXPERIMENT_STANDALONE)
  # experiment/ configured, built and tested on its own (with physics/), in a fresh build tree under this one.
  add_test(NAME experiment_standalone
           COMMAND ${CMAKE_COMMAND}
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}
                   -DBIN=${CMAKE_CURRENT_BINARY_DIR}/standalone
                   -DGENERATOR=${CMAKE_GENERATOR}
                   -DCXX=${CMAKE_CXX_COMPILER}
                   -DNAME=experiment
                   -DTEST_REGEX=^experiment_
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/../physics/tests/physics_standalone.cmake)
  set_tests_properties(experiment_standalone PROPERTIES TIMEOUT 600)
endif()
