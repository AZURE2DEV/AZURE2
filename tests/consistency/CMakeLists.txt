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

# The skill checker (tools/check_skills.py, docs/upgrade/RULES.md §3): its unit
# test, and the check itself over .claude/skills/ (user-level skills skipped so
# the result does not depend on the machine).
add_test(NAME check_skills_unit
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/check_skills_test.py)
add_test(NAME check_skills
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_SOURCE_DIR}/tools/check_skills.py --user-skills none)

# tools/classify.py (the comparison of two engines' outputs): a file that cannot be read is never identical.
add_test(NAME tools_classify
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/classify_test.py)

# tools/classify_quadrature.py (the gate of the shared mesh, docs/upgrade/energy_grid.md section 7.1): the per-file
# rule, largest and rms distance to the fine reference, and the fine reference taken per file from a list of tags.
add_test(NAME tools_classify_quadrature
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/classify_quadrature_test.py)

# tools/runall.sh and tools/run_one.sh (the whole-tree runs): a tag is never mixed with the results of an earlier run.
add_test(NAME tools_runall
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/runall_test.py)
