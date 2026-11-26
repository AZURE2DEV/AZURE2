# convert/: everything that knows the old project formats (docs/upgrade/project_format3.md section 8).
#
# azure2_convert_lines: the old line parsers, header-only (include/azure2/convert/*.h but Convert.h).  The engine
#   still reads format 2 through them until the format-3 reader (step 3b) replaces those sections; it includes them
#   from here, so that each old grammar has one parser.  3b removes the engine's use of this target.
cmake_minimum_required(VERSION 3.16...4.0)

add_library(azure2_convert_lines INTERFACE)
target_include_directories(azure2_convert_lines INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# azure2_convert: the converter old -> format 3 (Convert.h), used only by AZURE2 --convert ... --to 3 and the binding
#   (pyazr.convert_azr), and CheckFormat3Text, which AZURE2 --validate and pyazr.validate use until the format-3
#   reader (step 3b) holds it.  It links io/ (the document, the record reader and writer), experiment/ (the table and
#   the validator) and physics/ (element symbols), and reads the <levels> pair columns with io/'s LevelRow
#   (format 3 keeps that grammar).  It includes no engine header (repo_consistency check 9c).
# TEMPORARY (step 3d): the old sections as the engine's set-up specs (SpecsFromOldLines.h), the engine's one way to read
# an old project.  It uses the line parsers and experiment/'s SetupSpecs.h, nothing of the engine.
add_library(azure2_convert_specs STATIC src/SpecsFromOldLines.cpp)
target_include_directories(azure2_convert_specs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(azure2_convert_specs PUBLIC azure2_experiment azure2_convert_lines)
target_compile_features(azure2_convert_specs PUBLIC cxx_std_17)
set_target_properties(azure2_convert_specs PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_convert_specs PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_library(azure2_convert STATIC
  src/Build.cpp
  src/BuildApparatus.cpp
  src/BuildData.cpp
  src/BuildSections.cpp
  src/Convert.cpp
  src/Format2.cpp
  src/LineRepair.cpp
  src/OldEra.cpp
  src/OldProject.cpp
  src/PotentialRecords.cpp)
add_library(azure2::convert ALIAS azure2_convert)
target_include_directories(azure2_convert PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                                          PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(azure2_convert PUBLIC azure2_io azure2_experiment azure2_experiment_reader azure2_physics azure2_convert_lines azure2_convert_specs)
target_compile_features(azure2_convert PUBLIC cxx_std_17)
set_target_properties(azure2_convert PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_convert PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Tests: each construct of section 8.1 on a minimal old line, its records written out in the test
# (convert_constructs); the committed conversions of a fixed set of test projects, byte for byte, idempotence and the
# refusals (convert_fixtures).  (Sources named in full: repo_consistency checks that every tests/*_test.cpp is here.)
add_executable(convert_constructs_test tests/constructs_test.cpp)
target_link_libraries(convert_constructs_test PRIVATE azure2_convert)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(convert_constructs_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/constructs_scratch)
add_test(NAME convert_constructs COMMAND convert_constructs_test ${CMAKE_CURRENT_BINARY_DIR}/constructs_scratch)
# Variants of old lines (docs/upgrade/qt_gui.md section 19.5): the unambiguous ones converted with a note, the others
# refused with the line's number, the fields expected and the readings proposed (convert_line_repair).
add_executable(convert_line_repair_test tests/line_repair_test.cpp)
target_link_libraries(convert_line_repair_test PRIVATE azure2_convert)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(convert_line_repair_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/line_repair_scratch)
add_test(NAME convert_line_repair COMMAND convert_line_repair_test ${CMAKE_CURRENT_BINARY_DIR}/line_repair_scratch)
# The forms of the old project file through its history (tests/convert/fixtures/history/README.md): each era's form
# converted with its meaning (records equal to the same project in a later era's form), the others refused with
# their reason (convert_history).
add_executable(convert_history_test tests/history_test.cpp)
target_link_libraries(convert_history_test PRIVATE azure2_convert)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(convert_history_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME convert_history COMMAND convert_history_test ${CMAKE_CURRENT_SOURCE_DIR}/..)
# AZURE2 --convert under a decimal-comma locale gives the bytes it gives under "C" (skipped without such a locale).
add_test(NAME convert_cli_locale COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/../tests/convert/locale_check.sh
         $<TARGET_FILE:AZURE2> ${CMAKE_CURRENT_SOURCE_DIR}/..)
set_tests_properties(convert_cli_locale PROPERTIES SKIP_RETURN_CODE 77)
add_executable(convert_fixtures_test tests/fixtures_test.cpp)
target_link_libraries(convert_fixtures_test PRIVATE azure2_convert)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(convert_fixtures_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME convert_fixtures COMMAND convert_fixtures_test ${CMAKE_CURRENT_SOURCE_DIR}/..)
# The converter's first step, format 1 -> format 2 (ConvertToFormat2; docs/upgrade/project.md section 4), on every
# project file under tests/ and examples/ (the old ones are the originals in tests/convert/fixtures/).
add_executable(convert_format2_test tests/format2_test.cpp)
target_link_libraries(convert_format2_test PRIVATE azure2_convert)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(convert_format2_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME convert_format2 COMMAND convert_format2_test ${CMAKE_CURRENT_SOURCE_DIR}/../tests ${CMAKE_CURRENT_SOURCE_DIR}/../examples)
# Spec equivalence (project_format3.md step 3b): the specs of every old test project from its old lines equal those of
# its conversion through the format-3 reading path, member by member.
add_executable(convert_spec_equivalence_test tests/spec_equivalence_test.cpp)
target_link_libraries(convert_spec_equivalence_test PRIVATE azure2_convert)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(convert_spec_equivalence_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME convert_spec_equivalence COMMAND convert_spec_equivalence_test ${CMAKE_CURRENT_SOURCE_DIR}/..)
