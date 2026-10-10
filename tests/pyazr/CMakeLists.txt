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

# apply_fit needs the engine for its last section, and skips it when the
# module is not built -- so it is registered unconditionally like the above.
add_test(NAME pyazr_apply_fit
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/apply_fit_test.py)
set_tests_properties(pyazr_apply_fit PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Channel order within a level must not change a result -- and used to segfault.
add_test(NAME pyazr_channel_order
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/channel_order_test.py)
set_tests_properties(pyazr_channel_order PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Per-segment / per-dataset chi-squared, and AZURE2's own objective -- the
# penalty is cross-checked against the engine's Total-Norm-Chi-Squared.
add_test(NAME pyazr_objective
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/objective_test.py)
set_tests_properties(pyazr_objective PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Inactive <levels> and <segmentsData> lines must not shift pyazr's level
# numbering or its per-segment penalties off the engine's.
add_test(NAME pyazr_inactive_lines
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/inactive_lines_test.py)
set_tests_properties(pyazr_inactive_lines PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# write_output_files must produce the same files the binary does.
add_test(NAME pyazr_output_files
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/output_files_test.py)
set_tests_properties(pyazr_output_files PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The sqrt(E) energy-shift block on a segment line: read, written, and (with
# the engine) counted and penalised exactly as AZURE2 does.
add_test(NAME pyazr_sqrt_shift
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/sqrt_shift_test.py)
set_tests_properties(pyazr_sqrt_shift PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Park parametrization (--use-park): Brune equivalence of the cross sections,
# analytic gradient and residual Jacobian against finite differences (the
# adjoint differentiates Park's level matrix and its J-dependent terms), and
# the J > 0 penalty with its gradient.
add_test(NAME pyazr_park_gradient
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/park_gradient_test.py)
set_tests_properties(pyazr_park_gradient PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")

# Parameter files are read by name: old layouts, the #parametrization tag,
# write-back by name, foreign files refused.
add_test(NAME pyazr_sav_readers
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/sav_readers_test.py)
set_tests_properties(pyazr_sav_readers PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Two sessions in one process share no state: a data session's residuals are
# unchanged by an extrapolation (or second data) session opened beside it.
add_test(NAME pyazr_shared_state
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/shared_state_test.py)
set_tests_properties(pyazr_shared_state PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Tests that take minutes rather than seconds (measured 2026-10-10 on a CRC
# login node at their thread counts: park_gradient 41 min, shared_state 22,
# objective 7, inactive_lines 6.5, channel_order 6).  CI excludes them with
# `ctest -LE slow`; run the whole set (`ctest`, or `ctest -L slow` for these
# alone) before changing pyazr or the engine paths they exercise.
set_tests_properties(pyazr_channel_order pyazr_objective pyazr_inactive_lines
                     pyazr_park_gradient pyazr_shared_state
                     PROPERTIES LABELS slow)
