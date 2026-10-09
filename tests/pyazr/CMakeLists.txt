# The .azr round-trip test.  Pure Python and no compiled engine needed --
# AzrModel parses the file itself -- so it runs wherever python3 exists.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "python3 not found; skipping the pyazr round-trip test")
    return()
endif()

# Tests that compare pyazr with the CLI take the binary from AZURE2_BIN: the
# one this build made, not the newest build*/src/AZURE2 of the checkout.
add_test(NAME pyazr_roundtrip
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/roundtrip_test.py)

# apply_fit needs the engine for its last section, and skips it when the
# module is not built -- so it is registered unconditionally like the above.
add_test(NAME pyazr_apply_fit
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/apply_fit_test.py)
set_tests_properties(pyazr_apply_fit PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Channel order within a level must not change a result -- and used to segfault.
add_test(NAME pyazr_channel_order
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/channel_order_test.py)
set_tests_properties(pyazr_channel_order PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Per-segment / per-dataset chi-squared, and AZURE2's own objective -- the
# penalty is cross-checked against the engine's Total-Norm-Chi-Squared.
add_test(NAME pyazr_objective
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/objective_test.py)
set_tests_properties(pyazr_objective PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Inactive <levels> and <segmentsData> lines must not shift pyazr's level
# numbering or its per-segment penalties off the engine's.
add_test(NAME pyazr_inactive_lines
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/inactive_lines_test.py)
set_tests_properties(pyazr_inactive_lines PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# write_output_files must produce the same files the binary does.
add_test(NAME pyazr_output_files
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/output_files_test.py)
set_tests_properties(pyazr_output_files PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# A free THM norm is profiled, in pyazr as in the CLI: chi2, parameters,
# residuals, derivatives and output files on tests/7Li_p_a.
add_test(NAME pyazr_thm_profiled
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_profiled_test.py)
set_tests_properties(pyazr_thm_profiled PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# The THM sub-point grid follows a narrow level moved by a fit of the formal
# parameters (Brune off), not only under Brune.
add_test(NAME pyazr_thm_nonbrune_grid
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_nonbrune_grid_test.py)
set_tests_properties(pyazr_thm_nonbrune_grid PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# save_fit writes a fixed width's physical value as the fit had it (the fixed
# amplitude's width moves with the level's other widths), so it reloads as
# the same model.
add_test(NAME pyazr_fixed_width_roundtrip
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/fixed_width_roundtrip_test.py)
set_tests_properties(pyazr_fixed_width_roundtrip PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# add_extrapolation(frame="lab"|"cm"): the engine evaluates the grid asked for.
add_test(NAME pyazr_extrapolation_frame
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/extrapolation_frame_test.py)
set_tests_properties(pyazr_extrapolation_frame PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# THM experiments (<thm> experiment[...]): shared profiled norm and background
# in AzrModel, chi2, output files and derivatives, against the CLI.
add_test(NAME pyazr_thm_experiment
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_experiment_test.py)
set_tests_properties(pyazr_thm_experiment PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Coulomb line shape of the THM spectator (experiment[...] lineshape=on):
# AzrModel, CLI == pyazr, thm_lineshape against the formulas, model ratio.
add_test(NAME pyazr_thm_lineshape
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_lineshape_test.py)
set_tests_properties(pyazr_thm_lineshape PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Spectator-momentum window of a THM experiment (experiment[...] ps=...):
# AzrModel, CLI == pyazr, thm_vertex against the formulas, linearity in the nodes.
add_test(NAME pyazr_thm_spectator_window
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_spectator_window_test.py)
set_tests_properties(pyazr_thm_spectator_window PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Distortion factor R(E) of a THM experiment (experiment[...] distortion=...):
# AzrModel, CLI == pyazr, thm_distortion against mpmath, model ratio, a table.
add_test(NAME pyazr_thm_distortion
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_distortion_test.py)
set_tests_properties(pyazr_thm_distortion PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Fixed-angle THM observable (experiment[...] theta=...): AzrModel, CLI ==
# pyazr, 0-180 == angle-integrated / 4 pi, symmetry of an identical exit.
add_test(NAME pyazr_thm_fixed_angle
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_fixed_angle_test.py)
set_tests_properties(pyazr_thm_fixed_angle PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Distorted-wave entrance vertex (experiment[...] vertexModel=dw): AzrModel and
# its refusals, CLI == pyazr, thm_vertex (plane-wave limit), model ratio ==
# vertex ratio.
add_test(NAME pyazr_thm_dw_vertex
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_dw_vertex_test.py)
set_tests_properties(pyazr_thm_dw_vertex PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Spectator-direction window (experiment[...] spectatorAngles=...): AzrModel and
# its refusals, CLI == pyazr for R(E) and the DW vertex, thm_distortion and
# thm_vertex with a window (the reported vertex is the one the model uses).
add_test(NAME pyazr_thm_spectator_angles
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_spectator_angles_test.py)
set_tests_properties(pyazr_thm_spectator_angles PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# THM coherent background (cbackground=): AzrModel, CLI == session, Jacobian
# of the R-matrix and cbkg columns against differences, save_fit.
add_test(NAME pyazr_thm_coherent_background
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_coherent_background_test.py)
set_tests_properties(pyazr_thm_coherent_background PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Model averaging (pyazr/modelavg.py): Akaike/BIC/chi2/flat weights, priors
# and the chi2 rescale, stat/model/total spreads against closed forms, the
# averaged .azr writer.  Pure Python but for its last (engine) part.
add_test(NAME pyazr_model_average
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/model_average_test.py)
set_tests_properties(pyazr_model_average PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                     TIMEOUT 300)

# scripts/thm_model_average.py end to end on tests/18O_p_a_thm: dry run,
# 2 radii x pw/dw (dw refused and skipped), weights, averaged .azr.
add_test(NAME pyazr_thm_model_average
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_model_average_test.py)
set_tests_properties(pyazr_thm_model_average PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                     TIMEOUT 600)

# vertex=constant: the S_c(E_lowest) memo follows the channel radius across
# sessions in one process and after set_channel_radius on a live session.
add_test(NAME pyazr_thm_vertex_memo_sessions
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_vertex_memo_sessions_test.py)
set_tests_properties(pyazr_thm_vertex_memo_sessions PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                     TIMEOUT 300)

# save_fit writes the fitted direct-segment norms into the snapshot (the
# default) with explicit prior_centre rows, so it reproduces the fit's
# chi-squared and priors when run alone; or keeps the segment lines
# (norms="nominal") with the fitted values in the .sav.  AzrModel.
# prior_centres / set_prior_centres.
add_test(NAME pyazr_save_fit_norms
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/save_fit_norms_test.py)
set_tests_properties(pyazr_save_fit_norms PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# write_output_files returns the project's own output directory.
add_test(NAME pyazr_output_dir
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/output_dir_test.py)
set_tests_properties(pyazr_output_dir PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# thm_vertex flags (warns / raises) the energies its window does not reach.
add_test(NAME pyazr_thm_vertex_reach
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_vertex_reach_test.py)
set_tests_properties(pyazr_thm_vertex_reach PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# calculate_energies & co. run no observed-to-formal transformation, so an
# RWA vector no longer prints "Denominator less than zero".
add_test(NAME pyazr_geometry_no_transform
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/geometry_no_transform_test.py)
set_tests_properties(pyazr_geometry_no_transform PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# pyazr.tabulate drops the <thm> keys that name data segments (experiment
# lines, weight[k]) along with the segments, so a THM project tabulates.
add_test(NAME pyazr_tabulate_thm
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/tabulate_thm_test.py)
set_tests_properties(pyazr_tabulate_thm PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# THM under Park's parametrization: the same model as Brune's from the same
# file, and the THM Jacobian rows (central differences in Park's amplitudes)
# against differences of the residuals and the chain rule through J.
add_test(NAME pyazr_thm_park
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_park_test.py)
set_tests_properties(pyazr_thm_park PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# The sqrt(E) energy-shift block on a segment line: read, written, and (with
# the engine) counted and penalised exactly as AZURE2 does.
add_test(NAME pyazr_sqrt_shift
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/sqrt_shift_test.py)
set_tests_properties(pyazr_sqrt_shift PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Park parametrization (--use-park): Brune equivalence of the cross sections,
# analytic gradient and residual Jacobian against finite differences (the
# adjoint differentiates Park's level matrix and its J-dependent terms), and
# the J > 0 penalty with its gradient.
add_test(NAME pyazr_park_gradient
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/park_gradient_test.py)
set_tests_properties(pyazr_park_gradient PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Parameter files are read by name: old layouts, the #parametrization tag,
# write-back by name, foreign files refused.
add_test(NAME pyazr_sav_readers
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/sav_readers_test.py)
set_tests_properties(pyazr_sav_readers PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# A test that cannot run its engine part -- numpy, scipy or the compiled
# _azure2 module missing (pybind11 absent at configure time builds no module),
# or no AZURE2 binary -- exits 77 after whatever pure-Python checks it has, so
# ctest lists it as Skipped rather than Passed.  An exit of 0 there once hid
# that CI had never run a single CLI == pyazr comparison.  Set on every test
# in this directory, so a new one only has to follow the convention; the
# numpy-free roundtrip test never exits 77.  Keep this block LAST: it applies
# to the tests registered above it only (sqrt_shift, park_gradient and
# sav_readers, added after it, went without it until October 2026).
get_property(_pyazr_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${_pyazr_tests} PROPERTIES SKIP_RETURN_CODE 77)
