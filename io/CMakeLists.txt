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
  src/LevelGrouping.cpp)
add_library(azure2::io ALIAS azure2_io)
target_include_directories(azure2_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(azure2_io PUBLIC cxx_std_17)
set_target_properties(azure2_io PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_io PRIVATE -Wall -Wextra -Wpedantic)
endif()

# The lossless document on synthetic input; byte-identical round trips of every project file under tests/ and
# examples/ and of variants made from each (docs/upgrade/project.md section 3, P1).
# (Sources named in full: repo_consistency checks that every tests/*_test.cpp is named here.)
foreach(_src convert_test.cpp document_test.cpp level_grouping_test.cpp roundtrip_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(io_${_t}_test tests/${_src})
  target_link_libraries(io_${_t}_test PRIVATE azure2_io)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(io_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endforeach()
add_test(NAME io_document COMMAND io_document_test)
add_test(NAME io_level_grouping COMMAND io_level_grouping_test)
add_test(NAME io_roundtrip COMMAND io_roundtrip_test ${CMAKE_CURRENT_SOURCE_DIR}/../tests ${CMAKE_CURRENT_SOURCE_DIR}/../examples)
# The format-1 -> format-2 converter on every project file (docs/upgrade/project.md section 4).
add_test(NAME io_convert COMMAND io_convert_test ${CMAKE_CURRENT_SOURCE_DIR}/../tests ${CMAKE_CURRENT_SOURCE_DIR}/../examples)
