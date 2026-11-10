# io/: the project model (docs/upgrade/project.md).  Like physics/ it uses nothing else of AZURE2: no engine, GSL,
# Qt or OpenMP.  Public headers in include/azure2/io/, sources in src/, tests in tests/.
cmake_minimum_required(VERSION 3.16...4.0)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(azure2_io LANGUAGES CXX)
  enable_testing()
endif()

add_library(azure2_io STATIC
  src/Convert.cpp
  src/Document.cpp
  src/LevelGrouping.cpp
  src/LevelRows.cpp
  src/AbundanceLines.cpp
  src/RecordReader.cpp)
add_library(azure2::io ALIAS azure2_io)
target_include_directories(azure2_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(azure2_io PUBLIC cxx_std_17)
set_target_properties(azure2_io PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_io PRIVATE -Wall -Wextra -Wpedantic)
endif()

# The plain record structs of format 3 (include/azure2/io/Record.h, header-only), for experiment/: the table there
# validates records without linking io/ (docs/upgrade/project_format3.md section 7.3).  repo_consistency (check 9b)
# lets experiment/ include that one header of io/ and no other.
add_library(azure2_io_record INTERFACE)
target_include_directories(azure2_io_record INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(azure2_io_record INTERFACE cxx_std_17)

# The lossless document on synthetic input; byte-identical round trips of every project file under tests/ and
# examples/ and of variants made from each (docs/upgrade/project.md section 3, P1).
# (Sources named in full: repo_consistency checks that every tests/*_test.cpp is named here.)
foreach(_src abundance_lines_test.cpp convert_test.cpp document_test.cpp level_grouping_test.cpp level_rows_test.cpp
             record_reader_test.cpp roundtrip_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(io_${_t}_test tests/${_src})
  target_link_libraries(io_${_t}_test PRIVATE azure2_io)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(io_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endforeach()
add_test(NAME io_document COMMAND io_document_test)
add_test(NAME io_level_grouping COMMAND io_level_grouping_test)
# The <levels> and <abundances> grammars (LevelRows.h, AbundanceLines.h): synthetic cases, and every such section of the
# project files against the engine's readers before they moved to io/ (docs/upgrade/cnuc_split.md, step (c)).
add_test(NAME io_level_rows COMMAND io_level_rows_test ${CMAKE_CURRENT_SOURCE_DIR}/../tests ${CMAKE_CURRENT_SOURCE_DIR}/../examples)
add_test(NAME io_abundance_lines COMMAND io_abundance_lines_test ${CMAKE_CURRENT_SOURCE_DIR}/../tests ${CMAKE_CURRENT_SOURCE_DIR}/../examples)
# The generic record reader of format 3 (docs/upgrade/project_format3.md sections 3.1-3.3, 7.3): synthetic sections,
# each syntax code, and round trips of the record sections of the format-3 fixtures of experiment/tests/data/.
file(GLOB _format3_fixtures ${CMAKE_CURRENT_SOURCE_DIR}/../experiment/tests/data/*.azr)
add_test(NAME io_record_reader COMMAND io_record_reader_test ${_format3_fixtures})
add_test(NAME io_roundtrip COMMAND io_roundtrip_test ${CMAKE_CURRENT_SOURCE_DIR}/../tests ${CMAKE_CURRENT_SOURCE_DIR}/../examples)
# The format-1 -> format-2 converter on every project file (docs/upgrade/project.md section 4).
add_test(NAME io_convert COMMAND io_convert_test ${CMAKE_CURRENT_SOURCE_DIR}/../tests ${CMAKE_CURRENT_SOURCE_DIR}/../examples)
