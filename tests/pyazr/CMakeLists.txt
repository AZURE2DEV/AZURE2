# The .azr round-trip test.  Pure Python and no compiled engine needed --
# AzrModel parses the file itself -- so it runs wherever python3 exists.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "python3 not found; skipping the pyazr round-trip test")
    return()
endif()

add_test(NAME pyazr_roundtrip
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/roundtrip_test.py)
