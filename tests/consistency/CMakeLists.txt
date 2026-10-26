# Repository consistency: CLI flags vs --help and docs, test registration,
# check.sh modes, pyazr version and exports.  Reads the source tree only, so it
# needs nothing built and runs wherever python3 exists.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "python3 not found; skipping the repository consistency test")
    return()
endif()

add_test(NAME repo_consistency
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/repo_consistency_test.py)
