# support/: plumbing that every module may use and that uses nothing (CLAUDE.md, "Layout").  So far one piece: the
# environment switches, the one table of them and the one reader of the process's environment
# (include/azure2/support/Environment.h; docs/upgrade/evaluation.md section 8).  The stopping library and the engine
# link it; physics/ does not (it uses nothing).
add_library(azure2_support STATIC src/Environment.cpp)
add_library(azure2::support ALIAS azure2_support)
target_include_directories(azure2_support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(azure2_support PUBLIC cxx_std_17)
set_target_properties(azure2_support PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_support PRIVATE -Wall -Wextra -Wpedantic)
endif()

#   support_environment  the table (names unique, every kind and type known) and the reader (a record made from given
#                        values: set, empty, "0", numbers; read once per process)
add_executable(support_environment_test tests/environment_test.cpp)
target_link_libraries(support_environment_test PRIVATE azure2_support)
add_test(NAME support_environment COMMAND support_environment_test)
set_tests_properties(support_environment PROPERTIES ENVIRONMENT "AZURE_PROFILE=1;AZURE_GRAD_CHECK_H=2.5e-7")
