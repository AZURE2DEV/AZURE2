# experiment/: what was measured and what the model must predict (docs/upgrade/data_model.md section 2.5): the model
# expressions, from step S6 the lab -> c.m. frames, and the table of format 3 (docs/upgrade/project_format3.md
# section 7: record kinds, fields, quantities, units, the validator and the schema export).  It uses physics/
# (masses, constants, element symbols) and io/'s plain record structs (azure2/io/Record.h, header-only: the target
# azure2_io_record) and nothing else of AZURE2: no other io/ header, engine (src/), GSL, Qt or OpenMP.  The rule is
# enforced three times, as for physics/: the library links only azure2_physics and azure2_io_record; the ctest
# experiment_standalone configures and builds this directory (with physics/ and io/) as a project of its own; and
# repo_consistency refuses any include but azure2/experiment/, azure2/physics/, azure2/io/Record.h and the C++
# standard library.  The tests also link azure2_io: they read their fixtures with io's RecordReader, as step 3's
# reader will.
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
  # io/: the record structs (azure2_io_record) for the library, the reader (azure2_io) for the tests.
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../io ${CMAKE_CURRENT_BINARY_DIR}/io EXCLUDE_FROM_ALL)
endif()

add_library(azure2_experiment STATIC
  src/ColumnsGrid.cpp
  src/Diagnostics.cpp
  src/Build.cpp
  src/BuildEffects.cpp
  src/Expression.cpp
  src/Frames.cpp
  src/Json.cpp
  src/ModelDatum.cpp
  src/LevelTables.cpp
  src/Quantities.cpp
  src/ReactionParser.cpp
  src/Schema.cpp
  src/SchemaJson.cpp
  src/Units.cpp
  src/Validate.cpp
  src/ValueParsers.cpp
  src/Values.cpp)
add_library(azure2::experiment ALIAS azure2_experiment)
# The reading path of a format-3 project (ProjectReader.h): the table above with io/'s document, record and <levels>
# readers.
add_library(azure2_experiment_reader STATIC src/ProjectReader.cpp src/LevelTablesOfRows.cpp src/Describe.cpp)
target_link_libraries(azure2_experiment_reader PUBLIC azure2_experiment azure2_io)
set_target_properties(azure2_experiment_reader PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_experiment_reader PRIVATE -Wall -Wextra -Wpedantic)
endif()
target_include_directories(azure2_experiment PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(azure2_experiment PUBLIC azure2_physics azure2_io_record)
target_compile_features(azure2_experiment PUBLIC cxx_std_17)
set_target_properties(azure2_experiment PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_experiment PRIVATE -Wall -Wextra -Wpedantic)
endif()

# expression: each node against the engine's expression written out in the test, bit for bit, and backward against
# finite differences (docs/upgrade/data_model.md step S8).  frames: the lab -> c.m. functions against the engine's
# expressions before S6, bit for bit, and against closed forms and Lorentz boosts (step S6).  The format-3 table:
# the fixtures of tests/data/ (T2 of section 3.11, the user stories of section 9, the examples of sections 3.5-3.10)
# validate with no diagnostic, each diagnostic code is triggered alone, every quantity row's units and required
# fields are exercised, the value parsers on their grammars, and Schema::ToJson equals the committed
# tests/schema/schema_v3.json (check 21).
# (Sources named in full: repo_consistency checks that every tests/*_test.cpp is named here.)
foreach(_src expression_test.cpp frames_test.cpp schema_test.cpp diagnostics_test.cpp quantities_test.cpp
             parsers_test.cpp build_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(experiment_${_t}_test tests/${_src})
  target_link_libraries(experiment_${_t}_test PRIVATE azure2_experiment_reader)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(experiment_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME experiment_${_t} COMMAND experiment_${_t}_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/data
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/../tests/schema/schema_v3.json)
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
