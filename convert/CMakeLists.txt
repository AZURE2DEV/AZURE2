# convert/: everything that knows the old project formats (docs/upgrade/project_format3.md section 8).
#
# azure2_convert_lines: the old line parsers, header-only (include/azure2/convert/*.h but Convert.h).  The engine
#   still reads format 2 through them until the format-3 reader (step 3b) replaces those sections; it includes them
#   from here, so that each old grammar has one parser.  3b removes the engine's use of this target.
cmake_minimum_required(VERSION 3.16...4.0)

add_library(azure2_convert_lines INTERFACE)
target_include_directories(azure2_convert_lines INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
