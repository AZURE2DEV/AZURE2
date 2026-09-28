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

# A free THM norm is profiled, in pyazr as in the CLI: chi2, parameters,
# residuals, derivatives and output files on tests/7Li_p_a.
add_test(NAME pyazr_thm_profiled
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_profiled_test.py)
set_tests_properties(pyazr_thm_profiled PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The THM sub-point grid follows a narrow level moved by a fit of the formal
# parameters (Brune off), not only under Brune.
add_test(NAME pyazr_thm_nonbrune_grid
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_nonbrune_grid_test.py)
set_tests_properties(pyazr_thm_nonbrune_grid PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# save_fit writes a fixed width's physical value as the fit had it (the fixed
# amplitude's width moves with the level's other widths), so it reloads as
# the same model.
add_test(NAME pyazr_fixed_width_roundtrip
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/fixed_width_roundtrip_test.py)
set_tests_properties(pyazr_fixed_width_roundtrip PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")
