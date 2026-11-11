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

# calculate_chi2 (physical parameters) is the same objective as calculate_chi2_rwa:
# the normalizations in its vector are applied (they were not before 2026-10-26),
# the Park / abundance / scattering-length terms are included, a rejected vector is +inf.
add_test(NAME pyazr_physical_norm
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/physical_norm_test.py)
set_tests_properties(pyazr_physical_norm PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A channel radius set with azure2.set_channel_radius survives every later rebuild of the session
# (another pair's radius, a nuclear-potential change, a mode switch); before 2026-11-02 each rebuild
# re-read the .azr with only its own call's radius.  The wrapper checks the read-back radius.
add_test(NAME pyazr_radius_override
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/radius_override_test.py)
set_tests_properties(pyazr_radius_override PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A session parses its .azr once: calculate_rwa / objective do not re-read the file per segment per
# call (before 2026-11-03: 83 ms per call from NFS on the 3H+d model instead of 8 ms).
add_test(NAME pyazr_azr_parse_cache
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/azr_parse_cache_test.py)
set_tests_properties(pyazr_azr_parse_cache PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The API reports the model value of an E1-only / E2-only segment, as chi^2 uses it (ModelPass.h; HANDOFF).
add_test(NAME pyazr_api_model_value
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/api_model_value_test.py)
set_tests_properties(pyazr_api_model_value PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A rebuild (set_channel_radius, a potential change) keeps no Coulomb memo of the old radius, reads no
# unchanged input file and prints the setup Notes once per session (before 2026-11-04: ~96 kB per new
# radius, every file re-read, every Note reprinted; the 3H+d free-radius MCMC reached 239 GB).
add_test(NAME pyazr_rebuild_resources
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/rebuild_resources_test.py)
set_tests_properties(pyazr_rebuild_resources PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The hybrid Coulomb functions with a vanishing nuclear potential equal the plain ones: the
# derivative at the channel radius is O(dr^4) on a grid that has the radius as a node (before
# 2026-10-27 the chord slope of the bracketing cell put 4e-3 to 4e-2 into S; audit item 16).
add_test(NAME pyazr_hybrid_derivative
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/hybrid_derivative_test.py)
set_tests_properties(pyazr_hybrid_derivative PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Parameter.wigner_limit is AZURE2's 3 hbar^2 / (2 mu a^2); teichmann_wigner() returns it unchanged
# (it multiplied by 1.5 again before 2026-10-27, audit item 19) and lane_thomas_wigner() gives hbar^2/(mu a^2).
add_test(NAME pyazr_wigner_convention
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/wigner_convention_test.py)
set_tests_properties(pyazr_wigner_convention PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Inactive <levels> and <segmentsData> lines must not shift pyazr's level
# numbering or its per-segment penalties off the engine's.
add_test(NAME pyazr_inactive_lines
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/inactive_lines_test.py)
set_tests_properties(pyazr_inactive_lines PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# An external-capture final state must be a listed level: a used pair raises AZURE2Error, an unused one warns
# and gets no hidden level (tests/ec_final_state_listed models).
add_test(NAME pyazr_ec_final_state
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/ec_final_state_test.py)
set_tests_properties(pyazr_ec_final_state PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

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

# The first evaluation of a session uses the shift it is given: a first request of 0 against a nonzero shift on the
# line returned the chi2 of the line's shift when the segment had sub-points (17O+a Kunz profile, 2026-10-29).
add_test(NAME pyazr_shift_first_call
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/shift_first_call_test.py)
set_tests_properties(pyazr_shift_first_call PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Fitted target density (densityfit / densitygroup on a <targetInt> line): the depth window and no factor of N
# against the bare cross section, groups, penalties, the Jacobian column, shift + density, refused lines.
add_test(NAME pyazr_target_density
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/target_density_test.py)
set_tests_properties(pyazr_target_density PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Energy-shift functions (shiftfn): read, written, and (with the engine) listed,
# penalised and differenced exactly as AZURE2 does.
add_test(NAME pyazr_shift_fn
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/shift_fn_test.py)
set_tests_properties(pyazr_shift_fn PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Normalization groups (normgroup): read, written, and (with the engine) the
# follower's norm fixed, the penalty once and the shared column differenced.
add_test(NAME pyazr_norm_group
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/norm_group_test.py)
set_tests_properties(pyazr_norm_group PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Data covariance matrices ("covariance <file>"): the block read and written; with the
# engine, chi2 = sum of the whitened residuals squared, the gradient against central
# differences and 2 J^T w, and the relvector model frozen and refrozen.
add_test(NAME pyazr_data_covariance
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/data_covariance_test.py)
set_tests_properties(pyazr_data_covariance PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Park parametrization (--use-park): Brune equivalence of the cross sections,
# analytic gradient and residual Jacobian against finite differences (the
# adjoint differentiates Park's level matrix and its J-dependent terms), and
# the J > 0 penalty with its gradient.
add_test(NAME pyazr_park_gradient
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/park_gradient_test.py)
set_tests_properties(pyazr_park_gradient PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                     TIMEOUT 3600)

# Park's dS/dE at the level energies (bound, a few MeV, 33 MeV above threshold)
# against mpmath, with and without the Coulomb table, in one process.
add_test(NAME pyazr_park_dsde
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/park_dsde_test.py)
set_tests_properties(pyazr_park_dsde PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Parameter files are read by name: old layouts, the #parametrization tag,
# write-back by name, foreign files refused.
add_test(NAME pyazr_sav_readers
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/sav_readers_test.py)
set_tests_properties(pyazr_sav_readers PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# A parameter file's fixed level parameters reach the engine (load_sav), so
# pyazr evaluates the same model as the CLI given the same file.
add_test(NAME pyazr_sav_fixed_params
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/sav_fixed_params_test.py)
set_tests_properties(pyazr_sav_fixed_params PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

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

# Hard-sphere J^pi groups (includeAllJPi=1) through pyazr's default transform
# route: the level-less groups used to give the transform a 0 x 0 eigenproblem.
add_test(NAME pyazr_hard_sphere_transform
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/hard_sphere_transform_test.py)
set_tests_properties(pyazr_hard_sphere_transform PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Beta-decay feeding amplitudes keep their sign through the physical transform
# (GetTransformParams returned |B| until 2026-10-12), and a save_fit snapshot
# re-reads to the same chi-squared in the Brune, Park and standard
# parametrizations.
add_test(NAME pyazr_beta_sign
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/beta_sign_test.py)
set_tests_properties(pyazr_beta_sign PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A parameter vector passed as a list (or a non-float64 array) must work like
# an ndarray: the cast's temporary used to be freed with the GIL released.
add_test(NAME pyazr_list_args
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/list_args_test.py)
set_tests_properties(pyazr_list_args PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Several compound nuclei in one project: the <abundances> block, natural
# segments, the abundance fit parameters and merge_projects, against the engine.
add_test(NAME pyazr_natural_target
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/natural_target_test.py)
set_tests_properties(pyazr_natural_target PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# fetch_levels: IAEA LiveChart, falling back to the NNDC ENSDF API when LiveChart
# is unreachable or answers with a challenge page.  Offline (canned answers).
add_test(NAME pyazr_nds_levels_fallback
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/nds_levels_fallback_test.py)


# The parameter table, _all_rwa and the free vector follow the engine's layout
# in data and extrapolation mode (test segments change the parameter count).
add_test(NAME pyazr_extrap_layout
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/extrap_layout_test.py)
set_tests_properties(pyazr_extrap_layout PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# THM (HOES) derivatives: the adjoint against central differences
# (AZURE2_THM_FD=1) and Richardson differences, every branch of the adjoint.
add_test(NAME pyazr_thm_adjoint
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_adjoint_test.py)
set_tests_properties(pyazr_thm_adjoint PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# fetch_exfor / reference: the local EXFOR master files (IAEA-NDS/exfor_master) when the
# Web-API is unreachable or answers with a challenge page.  Offline: a synthetic entry in
# data/exfor_master and canned web failures.
add_test(NAME pyazr_nds_local_exfor
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/nds_local_exfor_test.py)

# Ratio-to-Rutherford EXFOR data with a lab angle column: to_azr converts the lab angle to the c.m.
# frame and divides by the Jacobian, so the engine sees ratio x sigma_R,cm(theta_cm) (audit item 8);
# the observable guess follows the angle column's frame.  Offline fixture subentries T0001004/5; the
# engine round trip runs when the extension is importable.
add_test(NAME pyazr_nds_rutherford
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/nds_rutherford_test.py)

# Angular-coefficient data segments (isDiff 9): add_data_segment(observable="angular-coefficients",
# absolute=...), the datasets reader, and the engine evaluating such a segment on the model's own
# Legendre coefficients (tests/angular_coefficients models).
add_test(NAME pyazr_angular_coefficients
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/angular_coefficients_test.py)
set_tests_properties(pyazr_angular_coefficients PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Relativistic kinematics: AzrModel.relativistic, and an engine session on the n + 4He hard-sphere model of
# tests/relativistic against sigma = 4 pi sin^2(ka)/k^2 with k from the invariant mass.
add_test(NAME pyazr_relativistic
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/relativistic_test.py)
set_tests_properties(pyazr_relativistic PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Scattering lengths: the <scatteringLengths> block helpers, azure2.scattering_lengths against the
# tests/scattering_length reference values, the constraint chi2, its gradient and Jacobian rows.
add_test(NAME pyazr_scattering_length
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/scattering_length_test.py)
set_tests_properties(pyazr_scattering_length PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Model data (<modelData>): pyazr.model_data's records against hand values on T2, the term in penalties(), its
# gradient against central differences, and the binding's resonance_strength / statistical_factor / branching_ratio.
add_test(NAME pyazr_model_data
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/model_data_test.py)
set_tests_properties(pyazr_model_data PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Ratio to Rutherford (data_kinds.md section 4.3): the chi^2 of ratio data against the same data multiplied by
# sigma_R, with and without an energy shift, the ratio of two folds, the binding's Coulomb cross sections and the
# gradient, against tests/ratio_to_rutherford/truth.py (mpmath; skipped without it).
add_test(NAME pyazr_ratio_to_rutherford
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/ratio_to_rutherford_test.py)
set_tests_properties(pyazr_ratio_to_rutherford PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Units of format 3 and the quantities that are scales (project_format3.md sections 3.3, 7.2): data files in other
# units to the bit or to 1e-12, S-factor data against an mpmath conversion (with a shift, and the gradient), yields per
# uC against per particle, pyazr.sfactor_to_cross_section (mpmath; skipped without it).
add_test(NAME pyazr_units_scales
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/units_scales_test.py)
set_tests_properties(pyazr_units_scales PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77)

# Electron screening: the AzrModel helpers for the screening block of a <targetInt> line, the
# screening_<n> parameter and its prior, and its analytic gradient and Jacobian column.
add_test(NAME pyazr_screening
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/screening_test.py)
set_tests_properties(pyazr_screening PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The external-capture node table keys on the relativistic switch: a session evaluating a model and then its
# relativistic=1 variant gives the fresh-process relativistic objective (before 2026-10-27 the nodes computed
# without the switch were reused; checklist t3-25).
add_test(NAME pyazr_ec_relativistic
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/ec_relativistic_test.py)
set_tests_properties(pyazr_ec_relativistic PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" LABELS "slow")

# The hybrid nuclear potential belongs to the session: a second session's <potential> block or set_potential
# used to reconfigure every live session through the process singleton (2026-10-27, audit item 30).
add_test(NAME pyazr_potential_sessions
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/potential_sessions_test.py)
set_tests_properties(pyazr_potential_sessions PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The shift function of a hybrid pair is continuous through threshold: below it the decaying solution of the radial
# equation with the potential (spin_orbit_potential.md section 2.4), against an independent integration.
add_test(NAME pyazr_hybrid_threshold
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/hybrid_threshold_test.py)
set_tests_properties(pyazr_hybrid_threshold PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The spin-orbit term of a hybrid pair through pyazr: channel functions per J against an independent integration, the
# API's term equal to the record's, the refusals, external capture per J (tests/spin_orbit_potential, ec_ratio.csv).
add_test(NAME pyazr_spin_orbit_potential
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/spin_orbit_potential_test.py)
set_tests_properties(pyazr_spin_orbit_potential PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The per-pair semantics of the hybrid potential on tests/hybrid_potential (its own folder's script): a pair's setting
# beats the default, a pair without one inherits it, the API and the file give the same chi-squared, the Coulomb memo
# sees a change of the potential; the project as committed against expected/chiSquared.out (about 20 s).
add_test(NAME pyazr_hybrid_per_pair
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/../hybrid_potential/check_per_pair.py)
set_tests_properties(pyazr_hybrid_per_pair PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The hybrid potential as a format-3 record: its Coulomb radius against an independent integration (DOP853, R_C =
# 5 fm), set_nuclear_potential(coulomb_radius=...) equal to the record, save_fit writing the record.
add_test(NAME pyazr_potential_record
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/potential_record_test.py)
set_tests_properties(pyazr_potential_record PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# chiSquared.out written through the API carries the Park, abundance and scattering-length totals
# (2026-10-27, audit item 30).
# The fit's objective through pyazr (docs/upgrade/objective.md step 4): objective() = objective_and_grad()[0] = the sum
# of the engine's terms, the gradient against finite differences, on models with norm, shift-function and
# scattering-length terms.
add_test(NAME pyazr_objective_terms
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/objective_terms_test.py)
set_tests_properties(pyazr_objective_terms PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# pyazr.schema(): the table of project format 3 equals tests/schema/schema_v3.json and AZURE2 --schema
# (docs/upgrade/project_format3.md section 6.1, check 21).
add_test(NAME pyazr_schema
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/schema_test.py)
set_tests_properties(pyazr_schema PROPERTIES ENVIRONMENT "AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# pyazr.extrapolation_grid and pyazr.refine_grid: the engine's grid forms and adaptive refinement; a session refuses
# an adaptive grid (docs/upgrade/project_format3.md section 3.8).
add_test(NAME pyazr_extrapolation_grid
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/extrapolation_grid_test.py)
set_tests_properties(pyazr_extrapolation_grid PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# AZURE2 --validate and pyazr.validate (the fixtures of experiment/tests/data/, a file with several errors, a format-2
# file refused) and pyazr.convert_azr(..., to=3) against the committed conversion (docs/upgrade/project_format3.md 6.2, 8).
add_test(NAME pyazr_validate
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/validate_test.py)
set_tests_properties(pyazr_validate PROPERTIES ENVIRONMENT "AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# pyazr.convert_azr: an old project in the updated format (docs/upgrade/project.md section 4) -- format 2, idempotent,
# the removed sections, the same chi-squared as the original.
add_test(NAME pyazr_convert_azr
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/convert_azr_test.py)
set_tests_properties(pyazr_convert_azr PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# AzrModel groups <levels> lines as the engine does (io::LevelGrouper; docs/upgrade/project.md P2): interleaved lines
# give the same levels, in AzrModel and in the engine.
add_test(NAME pyazr_level_grouping
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/level_grouping_test.py)
set_tests_properties(pyazr_level_grouping PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# AzrModel keeps a file's own layout (docs/upgrade/project.md section 5, step 4): unedited round trips of every
# project file and of CRLF / no-blank-line / interleaved variants; a field edit changes only its lines.
# One mass table (docs/upgrade/project.md section 5, step 5): the engine's THM nuclei carry the nuclear masses of
# pyazr/data/mass_1.mas20, the file fill_masses reads.
add_test(NAME pyazr_thm_nuclide_masses
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_nuclide_masses_test.py)
set_tests_properties(pyazr_thm_nuclide_masses PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

add_test(NAME pyazr_azr_layout
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/azr_layout_test.py)
set_tests_properties(pyazr_azr_layout PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

add_test(NAME pyazr_output_penalties
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/output_penalties_test.py)
set_tests_properties(pyazr_output_penalties PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The transform to physical parameters (level energies, external widths) leaves nothing in the Coulomb memo
# (2026-10-28, the 12C(a,g) MCMC memory growth).
add_test(NAME pyazr_coulomb_memo_transform
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/coulomb_memo_transform_test.py)
set_tests_properties(pyazr_coulomb_memo_transform PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# pyazr.masses reads the AME2020 table (pyazr/data/mass_1.mas20) exactly, its separation energies equal the
# engine's (m1 + m2 - M) c^2, and AzrModel.fill_masses rewrites only a pair's mass and S columns (2026-10-29).
add_test(NAME pyazr_masses
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/masses_test.py)

# The API calls behind pyazr calculate an evaluation's points in one parallel pass, as the CLI's chi-squared does,
# with results bit-identical to the serial path (2026-10-30; before, they ran one point at a time on one thread).
# Sets its own OMP_NUM_THREADS=4 in the processes it starts.
add_test(NAME pyazr_parallel_points
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/parallel_points_test.py)
set_tests_properties(pyazr_parallel_points PROPERTIES LABELS "slow")

# A session that reads the external-capture integrals back from intEC.dat equals the one that computed them
# (2026-10-30: the file had 6 significant digits; tests/13N chi2 moved by 0.62 on read-back).
add_test(NAME pyazr_ec_cache_precision
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/ec_cache_precision_test.py)
set_tests_properties(pyazr_ec_cache_precision PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# residual_jacobian(x, skip=cols): held columns are not finite-differenced and come back zero, everything else
# bit-identical (2026-10-30, rmfit polish passes the columns its bounds fix).
add_test(NAME pyazr_jacobian_skip
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/jacobian_skip_test.py)
set_tests_properties(pyazr_jacobian_skip PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# residual_jacobian differences shift and density columns forward only, a self-contained segment's shift on its own
# rows (2026-10-30: each side re-solved the whole Jacobian, adjoints included).
add_test(NAME pyazr_jacobian_shift_columns
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/jacobian_shift_columns_test.py)
set_tests_properties(pyazr_jacobian_shift_columns PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A level whose widths are all zero changes nothing, also at a point exactly at its energy (the level matrix is
# singular there; the A-matrix level route dropped the whole J-group, the R-matrix route gave 0/0).
add_test(NAME pyazr_zero_width_level
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/zero_width_level_test.py)
set_tests_properties(pyazr_zero_width_level PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The --ignore-externals mask (docs/upgrade/external_capture.md sections 3 and 7): defined on the input (1e-20 eV keeps
# its external part), formal = Brune capture with a masked level coupled to another, and the chi2 gradient equal to
# central differences in both bases (the formal one missed the projector's dependence, 1.5e-2, before 2026-11-08).
add_test(NAME pyazr_formal_capture_mask
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/formal_capture_mask_test.py)
set_tests_properties(pyazr_formal_capture_mask PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A width given as 0 is not fitted, in every basis (docs/upgrade/external_capture.md section 9): the formal basis frees
# the Brune basis's parameters (it freed rounding-valued ones), --no-ignore-externals keeps a zero photon width fixed,
# and apply_fit places a formal-basis vector in strict mode.
add_test(NAME pyazr_zero_width_fixed
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/zero_width_fixed_test.py)
set_tests_properties(pyazr_zero_width_fixed PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# pyazr.whittaker_wave: the engine's GSL-free Whittaker function (docs/upgrade/coulomb.md section 1.4) against
# mpmath, the scaled form outside the double range, broadcasting, ValueError for invalid arguments.
add_test(NAME pyazr_whittaker_wave
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/whittaker_wave_test.py)
set_tests_properties(pyazr_whittaker_wave PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The formal -> Brune iteration inside an evaluation (docs/upgrade/level_parameters.md section 9): chi2 of
# beam_profile_kernel in the formal basis is smooth in width_9_2 (it jumped 40.44 -> 42.93 when a level missed a 1e-12
# MeV tolerance by rounding and kept its formal values).  ~3 min.
add_test(NAME pyazr_formal_brune_iteration
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/formal_brune_iteration_test.py)
set_tests_properties(pyazr_formal_brune_iteration PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" LABELS "slow")

# The R-matrix route (<config> first line false) at a pole of the R-matrix: tests/13N's thick-target capture puts a
# sub-point exactly at a level's energy (NaN / 1.2 % off before 2026-11-01); identical_pp_res and natural_neon hit an
# exact pole with Brune (LU failure / 84 points 0).  ~4 min.
add_test(NAME pyazr_rmatrix_pole
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/rmatrix_pole_test.py)
set_tests_properties(pyazr_rmatrix_pole PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" LABELS slow)

# A failed point is NaN, also in SUM / RATIO segments, with a warning (before 2026-11-01: 0 or a missing term,
# silently).  Uses the AZURE_TEST_FAIL_E hook.
add_test(NAME pyazr_failed_point
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/failed_point_test.py)
set_tests_properties(pyazr_failed_point PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# use_brune=False on a project with no standard-R-matrix equivalent (tests/7Li_p_ay): AZURE2Error with the reason,
# not "z is not finite in log_Gamma" (before 2026-11-01).
add_test(NAME pyazr_no_brune_transform
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/no_brune_transform_test.py)
set_tests_properties(pyazr_no_brune_transform PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Composite segments and the route (2026-11-01): a zero RATIO denominator is NaN with a warning (was 0); component
# scalings keep every digit and one that is not a number refuses the project (was "%f" and a silent 1.0); the R-matrix
# route with ignore_externals=True runs the A-matrix with a warning (was the R-matrix route, 13N up to 4.5 % off).
add_test(NAME pyazr_composite_terms
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/composite_terms_test.py)
set_tests_properties(pyazr_composite_terms PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Free-gas Doppler broadening: add_target_effect(doppler=T) and the engine's Doppler step through the API
# (tests/doppler_fgm models).
add_test(NAME pyazr_doppler
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/doppler_test.py)
set_tests_properties(pyazr_doppler PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The fitted sample thickness of transmission / yield segments ("transmission n vary sigma", thicknessgroup):
# AzrModel helpers, type-10 records, penalties()["thickness"], chi2_and_grad and residual_jacobian columns
# (tests/transmission_yield models).
add_test(NAME pyazr_thickness
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thickness_test.py)
set_tests_properties(pyazr_thickness PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Parameter kinds (include/ParameterKind.h, docs/upgrade/parameter_kinds.md): the engine's table of kinds, and on
# projects holding every kind each parameter's name rebuilt from its record, the record against parameter_info, and
# pyazr's Parameter.kind / sample / thm_* from the records.
add_test(NAME pyazr_parameter_kinds
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/parameter_kinds_test.py)
set_tests_properties(pyazr_parameter_kinds PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The light particle's excitation energy, column 32 of a <levels> line (2026-11-05): AzrModel keeps the column
# (and adds or drops it only when e1 changes to or from 0), Pair carries light_excitation (and still reads an
# older engine's 19-field record), and the excitation moved between the partners gives the same cross sections
# (tests/light_excitation models).
add_test(NAME pyazr_light_excitation
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/light_excitation_test.py)
set_tests_properties(pyazr_light_excitation PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The prior parameter covariance ("priorCovariance=<file>" in <config>, 2026-11-05): pyazr.priors (write_prior /
# read_prior, read_covariance_dat, prior_from_fit), AzrModel.prior_covariance, and through the engine
# penalties()["prior_covariance"], objective(), the residual_jacobian rows L^-1 (p - mu) against central differences,
# the marginalization of a fixed listed parameter and a refused file (tests/norm_groups model).
add_test(NAME pyazr_prior_covariance
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/prior_covariance_test.py)
set_tests_properties(pyazr_prior_covariance PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Upgrade identity checks (docs/upgrade/ARCHITECTURE.md section 10; list and status in tests/upgrade/README.md).
# They use the test projects tests/upgrade/T1, T3 or small models built by tests/upgrade/lib, on temporary copies.
#   coulomb_identity: F, G, F', G' of the engine against mpmath, the Wronskian, P, S and the hard-sphere phase.
#   unitarity: one open channel (n + 16O): the optical-theorem total equals the elastic cross section.
#   optical_theorem: T3 n + 7Be total = elastic + (n,p) + (n,alpha alpha).
#   breit_wigner: an isolated narrow 12C(p,gamma) level integrates to the Breit-Wigner area.
add_test(NAME pyazr_upgrade_coulomb_identity
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/upgrade_coulomb_identity_test.py)
set_tests_properties(pyazr_upgrade_coulomb_identity PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")
add_test(NAME pyazr_upgrade_unitarity
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/upgrade_unitarity_test.py)
set_tests_properties(pyazr_upgrade_unitarity PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
add_test(NAME pyazr_upgrade_optical_theorem
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/upgrade_optical_theorem_test.py)
set_tests_properties(pyazr_upgrade_optical_theorem PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
add_test(NAME pyazr_upgrade_breit_wigner
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/upgrade_breit_wigner_test.py)
set_tests_properties(pyazr_upgrade_breit_wigner PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
#   parametrization_roundtrip: T1, T3 observables equal in the Brune, formal (--no-brune) and Park bases.
#   thread_independence: T1, T3 calculate / extrapolate / rate / fit at 1 and 4 threads.
#   gradient: analytic chi^2 gradient of T1, T3 = Richardson central differences.
add_test(NAME pyazr_upgrade_parametrization_roundtrip
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/upgrade_parametrization_roundtrip_test.py)
set_tests_properties(pyazr_upgrade_parametrization_roundtrip PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
add_test(NAME pyazr_upgrade_thread_independence
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/upgrade_thread_independence_test.py)
set_tests_properties(pyazr_upgrade_thread_independence PROPERTIES ENVIRONMENT "AZURE2_BIN=$<TARGET_FILE:AZURE2>" LABELS slow)
add_test(NAME pyazr_upgrade_gradient
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/upgrade_gradient_test.py)
set_tests_properties(pyazr_upgrade_gradient PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# pyazr.Project, the project as records (docs/upgrade/project_format3.md sections 4, 6): round trip of every converted
# project, one edit = one line, rename/remove, the engine's refusals, the user stories built from scratch against the
# committed projects' outputs (the CLI's), help() from the table, describe, session.
add_test(NAME pyazr_project
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/project_test.py)
set_tests_properties(pyazr_project PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# save_fit of a format-3 project with norms="fitted" (the fitted nuisances written into their records), a changed
# potential of a format-1/2 project refused with the conversion command, set_nuclear_potential's matching radius and step.
add_test(NAME pyazr_save_fit_format3
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/save_fit_format3_test.py)
set_tests_properties(pyazr_save_fit_format3 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The recipe "write a project from Python" (docs/source/user_guide/project_from_python.rst, the azure2-setup skill):
# pyazr/examples/project_from_python.py runs, and its chi^2 is the committed story's.
add_test(NAME pyazr_project_recipe
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/project_recipe_test.py)
set_tests_properties(pyazr_project_recipe PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")


# pyazr.predictions(path): the expression trees of a project, as AZURE2 --print-predictions prints them
# (docs/upgrade/data_model.md section 3, step S10).
add_test(NAME pyazr_predictions
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/predictions_test.py)
set_tests_properties(pyazr_predictions PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                     SKIP_RETURN_CODE 77)

# The gradient of an E1-only / E2-only SUM with a component (the component's full sigma, as the forward reads it)
# against central differences of the objective, by both component routes (docs/upgrade/data_model.md S11).
add_test(NAME pyazr_composite_gradient
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/composite_gradient_test.py)
set_tests_properties(pyazr_composite_gradient PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77
                     TIMEOUT 900)

# azure2.ec_integrals at the energy the engine uses for the final state (ALevel::GetECEnergy) in every basis: the
# formal basis gives the Brune basis's integrals on tests/13N (they were at the formal energy before 2026-11-10).
add_test(NAME pyazr_ec_integrals_formal
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/ec_integrals_formal_test.py)
set_tests_properties(pyazr_ec_integrals_formal PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")


# pyazr on a format-3 project without the old file class: the session's data sets from the engine's description, rename
# with the engine's parameter names, set keeping comments between repeats, records_from_csv.
add_test(NAME pyazr_session_records
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/session_records_test.py)
set_tests_properties(pyazr_session_records PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The nuclear part of a format-3 project in pyazr.Project (docs/upgrade/pyazr_on_records.md, "Edits of the nuclear
# part"): <levels> rows as levels and pairs, edits that change only the tokens they name, add / remove a level, the
# refusals, <config> entries, <thm> options and experiments, and save_fit's level values on 13N.
add_test(NAME pyazr_project_nuclear_part
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/project_nuclear_part_test.py)
set_tests_properties(pyazr_project_nuclear_part PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# pyazr.stopping: the binding of stopping/ (SRIM-2013 through LunaSRIM) against LunaSRIM's reference values
# (stopping/tests/data/srim_reference.csv), the same tolerances as the C++ tests.
add_test(NAME pyazr_stopping
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/stopping_test.py)
set_tests_properties(pyazr_stopping PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# A test that cannot run its engine part -- numpy, scipy or the compiled
# _azure2 module missing (pybind11 absent at configure time builds no module),
# or no AZURE2 binary -- exits 77 after whatever pure-Python checks it has, so
# ctest lists it as Skipped rather than Passed.  An exit of 0 there once hid
# that CI had never run a single CLI == pyazr comparison.  Set on every test
# in this directory, so a new one only has to follow the convention; the
# numpy-free roundtrip test never exits 77.  Keep this block LAST: it applies
# to the tests registered above it only (sqrt_shift, park_gradient and
# sav_readers went without it until October 2026; six more, the last
# pyazr_stopping, until November 2026: without PYAZR_ENGINE_DIR in the
# caller's environment they were skipped in a scratch build).
get_property(_pyazr_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${_pyazr_tests} PROPERTIES SKIP_RETURN_CODE 77)
# A build with -DPYAZR_OUTPUT_DIR outside pyazr/ (sanitizer and scratch builds)
# keeps _azure2 out of the checkout; PYAZR_ENGINE_DIR points the package at it,
# so these tests run the module this build made instead of skipping.
if(PYAZR_OUTPUT_DIR)
    set_property(TEST ${_pyazr_tests} APPEND PROPERTY ENVIRONMENT "PYAZR_ENGINE_DIR=${PYAZR_OUTPUT_DIR}")
endif()
