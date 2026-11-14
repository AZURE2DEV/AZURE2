# The self-tests of the Claude Code skills (.claude/skills/*/scripts/selftest.py; docs/upgrade/RULES.md section 3):
# each runs every script of its skill on a small project of tests/ in a temporary directory and prints PASS / FAIL
# lines.  One test per skill, named skill_<skill>; a selftest_slow.py beside it is skill_<skill>_slow with the
# label slow.  The environment is that of the pyazr tests: this build's AZURE2 and engine module, the checkout on
# PYTHONPATH (the skills import pyazr from it), one thread.  A self-test that cannot import the engine exits 77
# (Skipped).  AZURE2_COULOMB_TABLE comes from the caller's environment, as for the pyazr tests.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "python3 not found; skipping the skills' self-tests")
    return()
endif()

file(GLOB _skill_selftests CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/.claude/skills/*/scripts/selftest.py
                                            ${CMAKE_SOURCE_DIR}/.claude/skills/*/scripts/selftest_slow.py)
foreach(_selftest ${_skill_selftests})
    get_filename_component(_scripts ${_selftest} DIRECTORY)
    get_filename_component(_skill_dir ${_scripts} DIRECTORY)
    get_filename_component(_skill ${_skill_dir} NAME)
    get_filename_component(_kind ${_selftest} NAME_WE)
    set(_name skill_${_skill})
    if(_kind STREQUAL "selftest_slow")
        set(_name skill_${_skill}_slow)
    endif()
    add_test(NAME ${_name} COMMAND ${Python3_EXECUTABLE} ${_selftest})
    set(_env "OMP_NUM_THREADS=1;PYTHONPATH=${CMAKE_SOURCE_DIR};AZURE2_BIN=$<TARGET_FILE:AZURE2>")
    if(PYAZR_OUTPUT_DIR)
        list(APPEND _env "PYAZR_ENGINE_DIR=${PYAZR_OUTPUT_DIR}")
    endif()
    set_tests_properties(${_name} PROPERTIES ENVIRONMENT "${_env}" SKIP_RETURN_CODE 77)
    if(_kind STREQUAL "selftest_slow")
        set_tests_properties(${_name} PROPERTIES LABELS slow)
    endif()
endforeach()
