# The pyazr tests.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "python3 not found; skipping the pyazr tests")
    return()
endif()

# Tests that compare pyazr with the CLI take the binary from AZURE2_BIN: the
# one this build made, not the newest build*/src/AZURE2 of the checkout.

# save_fit's placement of level values (snapshot.place_level_values); needs the engine and skips (77) without it.
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

# the engine's chi-squared of a physical vector (the binding's calculate_chi2_physical) is that of the fit vector:
# the normalizations in its vector are applied (they were not before 2026-10-26),
# the Park / abundance / scattering-length terms are included, a rejected vector is +inf.
add_test(NAME pyazr_physical_norm
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/physical_norm_test.py)
set_tests_properties(pyazr_physical_norm PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A channel radius is a Project edit and a new session: the radii Project.set_pair writes are the session's, in
# either order of the edits, and a scan over radii in one process gives every radius its own model.
add_test(NAME pyazr_radius_override
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/radius_override_test.py)
set_tests_properties(pyazr_radius_override PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A channel's own radius is a channel_radius record added by Project.add and a new session: every l of a pair at one
# radius is the pair's radius edit; channel_functions with no radius is at the channel's own radius.
add_test(NAME pyazr_channel_radius_record
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/channel_radius_record_test.py)
set_tests_properties(pyazr_channel_radius_record PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Isospin levels across compound nuclei in a session: the shift and mix parameters and their priors, the objective and
# its rows, the gradient and the Jacobian through the ties.
add_test(NAME pyazr_isospin_levels
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/isospin_levels_test.py)
set_tests_properties(pyazr_isospin_levels PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A session parses its .azr once: evaluate().data[i].model / objective do not re-read the file per segment per
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

# Sessions of a radius scan in one process keep no Coulomb memo, read no unchanged input file again (the
# process's input-file cache) and each print their setup Notes once (a free-radius MCMC once reached 239 GB).
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

# Parameter.wigner_limit is AZURE2's 3 hbar^2 / (2 mu a^2), and pyazr.physics.wigner_limit gives the same number.
add_test(NAME pyazr_wigner_convention
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/wigner_convention_test.py)
set_tests_properties(pyazr_wigner_convention PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Session.wigner_widths = 2 P_l(E_r) gamma_W^2 against mpmath's Coulomb functions and the closed form of a neutral
# channel, T1, T2, T3 and 13N in the three parametrizations (docs/upgrade/pyazr3.md section 11); about 15 s.
add_test(NAME pyazr_wigner_widths
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/wigner_widths_test.py)
set_tests_properties(pyazr_wigner_widths PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# The Weisskopf estimate (physics/Widths.h through pyazr.weisskopf_width) against its closed form and against the
# numbers of the Python formula it replaced.
add_test(NAME pyazr_weisskopf
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/weisskopf_test.py)
set_tests_properties(pyazr_weisskopf PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# pyazr.endf reading: numbers, records, interpolation laws, MF=3, MF=6 LAW=2 and LAW=5 against hand values, curves at
# a data set's points, the column rules, the module's own format check, the ENDF/B-VIII.1 extract.
add_test(NAME pyazr_endf_read
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_read_test.py)
set_tests_properties(pyazr_endf_read PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# The model as a curve over a data set (Session.model_curve): at the data points equal to Session.evaluate on T1, T2
# and a ratio to Rutherford (unfolded to rounding, folded to the quadrature), between the points equal to the same
# points read as data, energy shifts, the grid, the refusals with their three parts (about 25 s).
add_test(NAME pyazr_model_curve
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/model_curve_test.py)
set_tests_properties(pyazr_model_curve PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# pyazr.endf's start model: the engine's two level doors against its own transformation (T1, T2, 13N), the
# ENDF/B-VIII.1 extract into a project reproducing its MF=3/MT=600, the refusals, write then read of T1 and T2.
add_test(NAME pyazr_endf_model
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_model_test.py)
set_tests_properties(pyazr_endf_model PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# pyazr.endf.write: T1, T2, 13N written and read back (own format check, MF=3 against the engine at the data and at
# the grid's midpoints, MF=2 back to the levels); the real checking codes print their findings when found.
# The quick set writes T1 and follows the narrow level; T2 and 13N (minutes: the elastic tables) are in the slow set.
add_test(NAME pyazr_endf_write
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_write_test.py T1 narrow)
set_tests_properties(pyazr_endf_write PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)
add_test(NAME pyazr_endf_write_slow
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_write_test.py T2 13N)
set_tests_properties(pyazr_endf_write_slow PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77
                     LABELS slow TIMEOUT 3000)
add_test(NAME pyazr_endf_write_tabulated
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_write_test.py T2 tabulated)
set_tests_properties(pyazr_endf_write_tabulated PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77
                     LABELS slow TIMEOUT 3000)

# MF=32 (LRF=7, LCOMP=1): the fit's covariance carried through the transformation to MF=2's formal parameters at T1's
# MIGRAD minimum (J against differences), written and read back, the refusals (about 45 s).
add_test(NAME pyazr_endf_covariance
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_covariance_test.py)
set_tests_properties(pyazr_endf_covariance PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# pyazr.endf.check_resonances on the committed 3He + 4He extract: the design's three findings (reconstruction and
# threshold, the 6Li mass against AME2020 and the Q values, the text against the records), and check(resonances=True).
add_test(NAME pyazr_endf_evaluation_check
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_evaluation_check_test.py)
set_tests_properties(pyazr_endf_evaluation_check PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# Two repairs of the ENDF export outside pyazr.endf: the Legendre frame of an extrapolation cut by a reshape (the
# arrays of the loop it replaces), Project.new without the parameter settings of removed records.
add_test(NAME pyazr_endf_export_support
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_export_support_test.py)
set_tests_properties(pyazr_endf_export_support PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# The engine's door of the nuclear amplitude expansion (LAW=5 LTP=1): the coefficients of T1, p + p and d + d give
# the engine's elastic cross section back by the manual's eqs. 6.9-6.14 (about 10 s).
add_test(NAME pyazr_endf_amplitudes
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_amplitudes_test.py)
set_tests_properties(pyazr_endf_amplitudes PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# pyazr.endf's checking codes against a stand-in (tests/pyazr/endf_standin/standin.py): finding, decks, reports read
# into findings, a crash as a finding, candidates and their difference, accept.  The real codes run too when
# AZURE2_ENDF_UTILITIES and AZURE2_PREPRO are set in the environment of ctest (skipped otherwise).
add_test(NAME pyazr_endf_check
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/endf_check_test.py)
set_tests_properties(pyazr_endf_check PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

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

# A capture final state above a particle threshold is refused: Project.validate() lists E_FINAL_STATE_UNBOUND and
# Session raises (T1 with its 3/2- final state above S_p); one just below the threshold opens.
add_test(NAME pyazr_ec_final_state_unbound
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/ec_final_state_unbound_test.py)
set_tests_properties(pyazr_ec_final_state_unbound PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# write_outputs must produce the same files the binary does.
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

# save_fit keeps a width stated as zero at zero when external capture gives its channel a physical width, so the
# snapshot has the fit's free parameters and no more.
add_test(NAME pyazr_save_fit_zero_width
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/save_fit_zero_width_test.py)
set_tests_properties(pyazr_save_fit_zero_width PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# an extrapolation grid in frame lab or cm: the engine evaluates the grid asked for.
add_test(NAME pyazr_extrapolation_frame
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/extrapolation_frame_test.py)
set_tests_properties(pyazr_extrapolation_frame PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# THM experiments (<thm> experiment[...]): shared profiled norm and background
# by Project edits, chi2, output files and derivatives, against the CLI.
add_test(NAME pyazr_thm_experiment
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_experiment_test.py)
set_tests_properties(pyazr_thm_experiment PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Coulomb line shape of the THM spectator (experiment[...] lineshape=on):
# Project edits, CLI == pyazr, thm_lineshape against the formulas, model ratio.
add_test(NAME pyazr_thm_lineshape
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_lineshape_test.py)
set_tests_properties(pyazr_thm_lineshape PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Spectator-momentum window of a THM experiment (experiment[...] ps=...):
# Project edits, CLI == pyazr, thm_vertex against the formulas, linearity in the nodes.
add_test(NAME pyazr_thm_spectator_window
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_spectator_window_test.py)
set_tests_properties(pyazr_thm_spectator_window PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Distortion factor R(E) of a THM experiment (experiment[...] distortion=...):
# Project edits, CLI == pyazr, thm_distortion against mpmath, model ratio, a table.
add_test(NAME pyazr_thm_distortion
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_distortion_test.py)
set_tests_properties(pyazr_thm_distortion PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Fixed-angle THM observable (experiment[...] theta=...): Project edits, CLI ==
# pyazr, 0-180 == angle-integrated / 4 pi, symmetry of an identical exit.
add_test(NAME pyazr_thm_fixed_angle
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_fixed_angle_test.py)
set_tests_properties(pyazr_thm_fixed_angle PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Distorted-wave entrance vertex (experiment[...] vertexModel=dw): Project edits and
# its refusals, CLI == pyazr, thm_vertex (plane-wave limit), model ratio ==
# vertex ratio.
add_test(NAME pyazr_thm_dw_vertex
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_dw_vertex_test.py)
set_tests_properties(pyazr_thm_dw_vertex PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Spectator-direction window (experiment[...] spectatorAngles=...): Project edits and
# its refusals, CLI == pyazr for R(E) and the DW vertex, thm_distortion and
# thm_vertex with a window (the reported vertex is the one the model uses).
add_test(NAME pyazr_thm_spectator_angles
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_spectator_angles_test.py)
set_tests_properties(pyazr_thm_spectator_angles PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# THM coherent background (cbackground=): Project edits, CLI == session, Jacobian
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
# sessions in one process, one after the other and alive at once.
add_test(NAME pyazr_thm_vertex_memo_sessions
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_vertex_memo_sessions_test.py)
set_tests_properties(pyazr_thm_vertex_memo_sessions PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                     TIMEOUT 300)

# save_fit writes the fitted direct-segment norms into the snapshot (the
# default) with explicit prior_centre rows, so it reproduces the fit's
# chi-squared and priors when run alone; or keeps the segment lines
# (norms="nominal") with the fitted values in the .sav.
# prior_centres / set_prior_centres.
add_test(NAME pyazr_save_fit_norms
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/save_fit_norms_test.py)
set_tests_properties(pyazr_save_fit_norms PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# write_outputs returns the project's own output directory.
add_test(NAME pyazr_output_dir
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/output_dir_test.py)
set_tests_properties(pyazr_output_dir PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# thm_vertex flags (warns / raises) the energies its window does not reach.
add_test(NAME pyazr_thm_vertex_reach
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_vertex_reach_test.py)
set_tests_properties(pyazr_thm_vertex_reach PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# evaluate().data[i].energies & co. run no observed-to-formal transformation, so an
# RWA vector no longer prints "Denominator less than zero".
add_test(NAME pyazr_geometry_no_transform
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/geometry_no_transform_test.py)
set_tests_properties(pyazr_geometry_no_transform PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# pyazr.tabulate drops the THM experiments and weights along with the data
# sets, so a THM project tabulates.
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

# chi^2 of a folded data set is smooth in its energy shift: the kernel moves on the fixed nodes of the shared mesh,
# clipped at its moved support (include/ShiftedFold.h); 21 shifts across 0.2 keV in one session.
add_test(NAME pyazr_shift_fold_smooth
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/shift_fold_smooth_test.py)
set_tests_properties(pyazr_shift_fold_smooth PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# A thick target with a positive energy shift on the fixed nodes: the translated stopping is never evaluated below the
# nominal support (no "unsupported imaginary number" from a power-law stopping), and the yields equal those of the
# same shift planned for exactly (include/ShiftedFold.h).
add_test(NAME pyazr_shift_fold_floor
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/shift_fold_floor_test.py)
set_tests_properties(pyazr_shift_fold_floor PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# A folded capture data set with external capture follows its shift beyond the reach: the re-planned sites take their
# capture amplitudes from the tables, no refusal mid-fit.
add_test(NAME pyazr_shift_fold_ec
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/shift_fold_ec_test.py)
set_tests_properties(pyazr_shift_fold_ec PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

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
# The same with the node tables checked at every hand-over and point pass (AZURE_NODE_CHECK=1): the Jacobian's local
# shift column moves the points of its segments and must look their slots up again (FillShiftColumns), as the
# hand-over does; with the slots left on the old nodes objective().residuals_and_jacobian raised "48 differ".
add_test(NAME pyazr_shift_fn_node_check
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/shift_fn_test.py)
set_tests_properties(pyazr_shift_fn_node_check PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE_NODE_CHECK=1")

# Normalization groups (normgroup): read, written, and (with the engine) the
# follower's norm fixed, the penalty once and the shared column differenced.
add_test(NAME pyazr_norm_group
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/norm_group_test.py)
set_tests_properties(pyazr_norm_group PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Data covariance matrices ("covariance <file>"): the block read and written; with the
# engine, chi2 = sum of the whitened residuals squared, the gradient against central
# differences and 2 J^T w, and the relvector model frozen and refrozen.
# The rate libraries' formats (pyazr/rate_formats.py): a REACLIB fit of known sets, the files written and read back,
# STARLIB with a band on T1, the refusals.
add_test(NAME pyazr_rate_formats
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/rate_formats_test.py)
set_tests_properties(pyazr_rate_formats PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77)

# The data covariance builder (pyazr/covariance.py): parts against matrices worked by hand, its refusals, the file it
# writes, and the engine's chi-squared with the file equal to v^T C^-1 v.
add_test(NAME pyazr_covariance_builder
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/covariance_builder_test.py)
set_tests_properties(pyazr_covariance_builder PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77)

# The reaction rate through pyazr (pyazr/rate.py): the program's rates digit for digit on a grid, a list and with a
# parameter file; exit -1; the bands against direct loops; the Gamow window against its closed form.
add_test(NAME pyazr_reaction_rate
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/reaction_rate_test.py)
set_tests_properties(pyazr_reaction_rate PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                     SKIP_RETURN_CODE 77)

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

# Two levels of one J^pi at one energy (physics/include/azure2/physics/LevelPairShift.h): chi-squared continuous
# across the gap in the Brune and Park bases, the analytic gradient against central differences at the coincidence
# and beside it, no jump along a path through it.
add_test(NAME pyazr_coincident_levels
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/coincident_levels_test.py)
set_tests_properties(pyazr_coincident_levels PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77
                     TIMEOUT 1800)

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
set_tests_properties(pyazr_sav_fixed_params PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Two sessions in one process share no state: a data session's residuals are
# unchanged by an extrapolation (or second data) session opened beside it.
add_test(NAME pyazr_shared_state
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/shared_state_test.py)
set_tests_properties(pyazr_shared_state PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# No hidden state (docs/upgrade/evaluation.md): projects A, B (and others) in one process in the orders A B A and
# B A B, at 1 and 3 threads, each evaluation byte-identical to the project alone in a fresh process.  About 1.5 min;
# its worker processes use up to 3 threads.
add_test(NAME pyazr_session_isolation
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/session_isolation_test.py)
set_tests_properties(pyazr_session_isolation PROPERTIES PROCESSORS 3 SKIP_RETURN_CODE 77)
# Two sessions evaluating at the same time on two Python threads (and two Coulomb table files in one process); the
# second table is made with the table tool of this build.
add_test(NAME pyazr_concurrent_sessions
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/concurrent_sessions_test.py)
set_tests_properties(pyazr_concurrent_sessions PROPERTIES PROCESSORS 2 SKIP_RETURN_CODE 77 LABELS slow
                     ENVIRONMENT "AZURE2_COULOMB_TABLE_TOOL=$<TARGET_FILE:azure2-coulomb-table>")

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

# Several compound nuclei in one project: a target with isotopes, the abundance fit parameters and merge_projects
# (the engine's merge of the <levels> rows), against the engine.
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
# Its own limit: 1081 s on the Windows CI machine, past ctest's 2400 s on the macOS one (76875f72).
set_tests_properties(pyazr_thm_adjoint PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 5400)

# fetch_exfor / reference: the local EXFOR master files (IAEA-NDS/exfor_master) when the
# Web-API is unreachable or answers with a challenge page.  Offline: a synthetic entry in
# data/exfor_master and canned web failures.
add_test(NAME pyazr_nds_local_exfor
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/nds_local_exfor_test.py)

# EXFOR Web-API requests send the user-agent the IAEA Nuclear Data Section gave AZURE2 (Exfor-Azure2/1.0);
# LiveChart, the NNDC ENSDF API and Crossref keep pyazr's own.  Offline: urlopen is replaced by a recorder.
add_test(NAME pyazr_nds_exfor_user_agent
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/nds_exfor_user_agent_test.py)

# Ratio-to-Rutherford EXFOR data with a lab angle column: to_azr converts the lab angle to the c.m.
# frame and divides by the Jacobian, so the engine sees ratio x sigma_R,cm(theta_cm) (audit item 8);
# the observable guess follows the angle column's frame.  Offline fixture subentries T0001004/5; the
# engine round trip runs when the extension is importable.
add_test(NAME pyazr_nds_rutherford
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/nds_rutherford_test.py)

# Angular-coefficient data sets (quantity legendre, expansion a_over_a0 / absolute): the record, and the engine
# evaluating such a data set on the model's own
# Legendre coefficients (tests/angular_coefficients models).
add_test(NAME pyazr_angular_coefficients
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/angular_coefficients_test.py)
set_tests_properties(pyazr_angular_coefficients PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Relativistic kinematics: Project.config()/set_config(relativistic=), and an engine session on the n + 4He hard-sphere model of
# tests/relativistic against sigma = 4 pi sin^2(ka)/k^2 with k from the invariant mass.
add_test(NAME pyazr_relativistic
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/relativistic_test.py)
set_tests_properties(pyazr_relativistic PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Scattering lengths: the <scatteringLengths> block helpers, Session.scattering_lengths against the
# tests/scattering_length reference values, the constraint chi2, its gradient and Jacobian rows.
add_test(NAME pyazr_scattering_length
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/scattering_length_test.py)
set_tests_properties(pyazr_scattering_length PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Model data (<modelData>): pyazr.model_data's records against hand values on T2, the term in evaluate().terms, its
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

# Electron screening: the target record's screening field (Project), the
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

# A pair's potential through the session (potential_background.md, review 2026-11-20): its `wave` fields, the channel
# functions J-aware and dressed together, the potential's box states and its bound state against scipy.
add_test(NAME pyazr_potential_api
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/potential_api_test.py)
set_tests_properties(pyazr_potential_api PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77)

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
# 5 fm), a Project edit of coulomb_radius equal to the record, save_fit keeping the record.
add_test(NAME pyazr_potential_record
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/potential_record_test.py)
set_tests_properties(pyazr_potential_record PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# chiSquared.out written through the API carries the Park, abundance and scattering-length totals
# (2026-10-27, audit item 30).
# The fit's objective through pyazr (docs/upgrade/objective.md step 4): objective() = objective().value_and_grad[0] = the sum
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

# Project.levels() groups <levels> rows as the engine does (io::LevelGrouper): interleaved rows give the same levels,
# in Project and in the engine.
add_test(NAME pyazr_level_grouping
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/level_grouping_test.py)
set_tests_properties(pyazr_level_grouping PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Project keeps a file's own layout (the lossless document): unedited round trips of every
# project file and of CRLF / no-blank-line / interleaved variants; a field edit changes only its lines.
# One mass table (docs/upgrade/project.md section 5, step 5): the engine's THM nuclei carry the nuclear masses of
# pyazr/data/mass_1.mas20, the file pyazr.sources.masses reads.
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

# pyazr.sources.masses reads the AME2020 table (pyazr/data/mass_1.mas20) exactly, its separation energies equal the
# engine's (m1 + m2 - M) c^2, and masses.fill_masses (Project.set_pair) rewrites only a pair's masses and S (2026-10-29).
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

# objective().residuals_and_jacobian(x, skip=cols): held columns are not finite-differenced and come back zero, everything else
# bit-identical (2026-10-30, rmfit polish passes the columns its bounds fix).
add_test(NAME pyazr_jacobian_skip
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/jacobian_skip_test.py)
set_tests_properties(pyazr_jacobian_skip PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# objective().residuals_and_jacobian differences shift and density columns forward only, a self-contained segment's shift on its own
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

# pyazr.physics.collision_matrix: the level matrix and the collision matrix of one J^pi block at one energy from
# arrays, without a session (physics/CollisionMatrix.h through api/src/pybind_collision_matrix.cpp), against the
# level matrix written in numpy (formal, Brune's and Park's parameters, photon channels), unitarity, the one-level
# phase shift and Breit-Wigner form, the three parametrizations of one block giving one U, and the errors.
add_test(NAME pyazr_collision_matrix
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/collision_matrix_test.py)
set_tests_properties(pyazr_collision_matrix PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# The same function against a session's own numbers on T1, T2 and T3 in the three parametrizations: the arrays of
# every block read off the session, phase shifts and cross sections formed from U and T in numpy, compared with
# Session.evaluate (the capture cross sections on a copy of T1 without external capture).
add_test(NAME pyazr_collision_matrix_session
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/collision_matrix_session_test.py)
set_tests_properties(pyazr_collision_matrix_session PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# pyazr.physics.optical_scattering: elastic scattering in an optical potential (physics/transfer through
# api/src/pybind_optical.cpp) against FRESCO's reference files of physics/tests/data/fresco, closed forms, the errors.
add_test(NAME pyazr_optical_scattering
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/optical_scattering_test.py)
set_tests_properties(pyazr_optical_scattering PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# pyazr.physics: every row of the binding's physics table (api/physics_table/src/PhysicsTable.cpp) against mpmath, the engine's
# own use on 13N and the Python formula it replaces (the largest difference per formula is printed); broadcasting,
# errors="mask", DomainError (docs/upgrade/pyazr3.md section 2.7).  tests/consistency check 41 asks for a row_test per row.
add_test(NAME pyazr_physics_table
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/physics_table_test.py)
set_tests_properties(pyazr_physics_table PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

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

# Free-gas Doppler broadening: the engine's Doppler step through the API and the fitted temperature (<target>.doppler)
# (tests/doppler_fgm models, format 3).
add_test(NAME pyazr_doppler
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/doppler_test.py)
set_tests_properties(pyazr_doppler PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The fitted sample thickness of transmission / yield data sets (a sample target's fitted areal_density):
# type-10 records, terms["thickness priors"], objective().value_and_grad and objective().residuals_and_jacobian columns
# (tests/transmission_yield models, format 3).
add_test(NAME pyazr_thickness
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thickness_test.py)
set_tests_properties(pyazr_thickness PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Parameter kinds (inference/include/azure2/inference/ParameterKind.h, docs/upgrade/parameter_kinds.md): the engine's
# table of kinds, and on
# projects holding every kind each parameter's name rebuilt from its record, the records' ordinals and keys, and
# pyazr's Parameter.kind / sample / thm_* from the records.
add_test(NAME pyazr_parameter_kinds
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/parameter_kinds_test.py)
set_tests_properties(pyazr_parameter_kinds PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The parameter table and the pairs from the engine's records (api/src/pybind_records.cpp, pyazr.Parameters): every
# kind described by an owner field of the record table; on T2, 13N, the gamma_decay models and a THM project with a
# coherent background the names and start values against the program's param.par, the fixed flags against the
# project, x0 / free / of_kind / of_level / put / take / read / write, the pair records against the project's pairs.
add_test(NAME pyazr_parameter_records
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/parameter_records_test.py)
set_tests_properties(pyazr_parameter_records PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# A level's radiative decay as Gamma_gamma, branchings and mixing ratios (the gamma_decay record,
# docs/upgrade/gamma_parametrisation.md; tests/gamma_decay/models): the record built in Python, Gauss-Newton fits in
# both parametrisations (the same minimum, the branchings' errors equal the propagated ones) and an agent's
# read-change-refit.  Needs scipy (skipped without it).
add_test(NAME pyazr_gamma_decay
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/gamma_decay_test.py)
set_tests_properties(pyazr_gamma_decay PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77)

# The light particle's excitation energy, column 32 of a <levels> line (2026-11-05): Project keeps the column
# (Project.set_pair(light_excitation=) writes it, a row without it gets it only for a value other than 0), Pair carries light_excitation (and still reads an
# older engine's 19-field record), and the excitation moved between the partners gives the same cross sections
# (tests/light_excitation models).
add_test(NAME pyazr_light_excitation
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/light_excitation_test.py)
set_tests_properties(pyazr_light_excitation PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The prior parameter covariance ("priorCovariance=<file>" in <config>, 2026-11-05): pyazr.priors (write_prior /
# read_prior, read_covariance_dat, prior_from_fit), Project.config()/set_config(priorCovariance=), and through the engine
# terms["prior covariance"], objective(), the objective().residuals_and_jacobian rows L^-1 (p - mu) against central differences,
# the marginalization of a fixed listed parameter and a refused file (tests/norm_groups model).
add_test(NAME pyazr_prior_covariance
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/prior_covariance_test.py)
set_tests_properties(pyazr_prior_covariance PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The prior builder (pyazr.priors.build): the reader of the prior file against the engine's verdicts on a table of
# small files, a built prior's term through a session = d^T C^-1 d, the part fit = prior_from_fit, the part chain = the
# draws' mean and covariance, every refusal and note (tests/norm_groups model).
add_test(NAME pyazr_prior_builder
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/prior_builder_test.py)
set_tests_properties(pyazr_prior_builder PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77)

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
# potential of a format-1/2 project refused with the conversion command, a record's matching radius and step.
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

# Session.ec_integrals at the energy the engine uses for the final state (ALevel::GetECEnergy) in every basis: the
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

# The THM columns of <levels> as Project edits: set_pair(binding_energy=) and set_channel(width_is_amplitude=) give
# the chi-squared of the same edit made by hand on two THM projects; the refusals.
add_test(NAME pyazr_thm_pair_edits
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_pair_edits_test.py)
set_tests_properties(pyazr_thm_pair_edits PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Pairs added and removed (Project.add_pair, remove_pair) and the channels a J^pi allows (allowed_channels, the Qt
# rule): T1's nuclear part built from nothing gives T1's chi-squared; add then remove gives the text back; refusals.
add_test(NAME pyazr_pair_edits
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/pair_edits_test.py)
set_tests_properties(pyazr_pair_edits PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The keys of <thm> as the engine's table (src/ThmKeys.cpp, pyazr.kinds()["thm"]): the refusals name the table's
# keys, the web page's table agrees, every key reads what the table offers, the reference page is the table's.
add_test(NAME pyazr_thm_keys
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/thm_keys_test.py)
set_tests_properties(pyazr_thm_keys PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The session's diagnostic calls described as rows (pyazr.session_functions): every key of the answers of
# channel_functions and the THM calls on four projects is a path of its row, and every path is in an answer.
add_test(NAME pyazr_session_functions
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/session_functions_test.py)
set_tests_properties(pyazr_session_functions PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Small defects, each with the check that failed before: save_fit to another directory reopens, `python -m pyazr
# --help` without the engine, `pyazr ui FILE --json` keeps serving, ChannelWidth.jpi, a pair record's is_entrance.
add_test(NAME pyazr_small_defects
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/small_defects_test.py)
set_tests_properties(pyazr_small_defects PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# pyazr.stopping: the binding of stopping/ (SRIM-2013 through LunaSRIM) against LunaSRIM's reference values
# (stopping/tests/data/srim_reference.csv), the same tolerances as the C++ tests.
add_test(NAME pyazr_stopping
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/stopping_test.py)
set_tests_properties(pyazr_stopping PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# pyazr's errors (pyazr/errors.py): a table of wrong calls on tests/13N, each raising its class (a pyazr.Error and the
# builtin it raised before) with a message that says what was wrong, what is allowed and which call helps; the
# engine's notice of an EC amplitude file that cannot be written is a pyazr.EngineWarning.
add_test(NAME pyazr_errors
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/errors_test.py)
set_tests_properties(pyazr_errors PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# pyazr.spin_operator (pyazr/spin.py): the operator of a spin observable for its name in the literature; A_ij by the
# spins or refused; every operator the table gives validates in tests/spin_observables_p3He and _dd.
add_test(NAME pyazr_spin_names
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/spin_names_test.py)
set_tests_properties(pyazr_spin_names PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" SKIP_RETURN_CODE 77)

# The package describes itself from the engine's tables: pyazr.describe(), pyazr.kinds() (record, data and parameter
# kinds against AZURE2 --schema / --describe, options against the binding's RuntimeOptions and AZURE2 --help),
# Session.describe() on tests/upgrade/T1, and the session of the package docstring runs.
add_test(NAME pyazr_describe
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/describe_test.py)
set_tests_properties(pyazr_describe PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# The nuisance priors belong to the session's model: a session of a project edited to another radius, made after
# another session evaluated its objective, gives the nuisance term 1 at the start, and the same numbers twice.
add_test(NAME pyazr_rebuild_nuisance
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/rebuild_nuisance_test.py)
set_tests_properties(pyazr_rebuild_nuisance PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# The web UI's server (pyazr/ui, docs/upgrade/web_ui.md): pytest files in tests/pyazr/ui, each run as a script that
# exits 77 (skipped) with its reason when fastapi, httpx, uvicorn, pytest or the engine module is missing (the last
# block gives them the skip code and PYAZR_ENGINE_DIR).  The injected errors that showed each rule's test failing are
# listed in tests/pyazr/ui/README.md.
# static_test: the page's files as served (token cookie, CSP, every loaded file tracked, plotly's hash, a widget per
# value type); node_test: the page's modules under node's own test runner (skipped where node is absent);
# results_test: the Results payload against the program's AZUREOut_* and chiSquared.out (T1, T2, T3);
# engine_jobs_test: the engine's fitters and MINOS as jobs against Session.fit_with_engine (T3), a stop on T2.
# sampler_bands_test: emcee and the bands as jobs on a two-parameter copy of T3 against the Laplace approximation
# (about 30 s).  thm_routes_test: the THM tab's routes on T2 and a THM copy of 18O_p_a_thm against the session's calls
# and the program, the refusals on T1.  nucleus_routes_test: the Nucleus tab's pair, binding-energy, amplitude-flag
# and new-level edits against the same Project calls (T1, T2).  utilities_test: the covariance builder's routes (the
# file and the field as one history entry) and the rate job against pyazr.covariance and pyazr.rate (T1).
# prior_builder_test: the prior builder's routes against pyazr.priors.build, the file and the key as one entry (T1).
# curves_test: the model curve route and an ENDF file on the curve's grid against pyazr.model_curve (T1, T2).
# endf_project_test: a new project from an ENDF-6 file and the project written as one, against pyazr.endf (T1).
# first_project_test: T1 made from nothing through the page's routes (an empty project, the pairs by name, the levels,
# the data sets, the calculation equal to T1's, a fit, Save; about 60 s).
foreach(_ui_test ui/api_test.py ui/jobs_test.py ui/server_test.py ui/routes_consistency_test.py ui/static_test.py
        ui/node_test.py ui/results_test.py ui/engine_jobs_test.py ui/sampler_bands_test.py ui/thm_routes_test.py
        ui/nucleus_routes_test.py ui/utilities_test.py ui/endf_test.py
        ui/prior_builder_test.py ui/curves_test.py ui/endf_project_test.py ui/first_project_test.py)
    string(REGEX REPLACE "ui/(.*)_test.py" "pyazr_ui_\\1" _name ${_ui_test})
    add_test(NAME ${_name} COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${_ui_test})
    set_tests_properties(${_name} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
endforeach()
# The browser smoke test of the page: a server on T1 and T2, headless Chrome runs the page's self-test and takes one
# screenshot per tab (tests/pyazr/ui/smoke_test.py; about 30 s; skipped where no Chrome or Chromium is found).
# Its own label: ctest -L browser.
add_test(NAME pyazr_ui_smoke COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ui/smoke_test.py)
set_tests_properties(pyazr_ui_smoke PROPERTIES LABELS browser
                     ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
# The navigation walked by clicks in headless Chrome over its DevTools protocol (tests/pyazr/ui/navigation_test.py:
# the tab row in the old program's order, each main tab's table, a level chosen into the address, Reload and Back,
# Calculate, the validator's list, a fit followed and stopped from the status line; about 60 s; skipped without
# Chrome or the package websockets).
add_test(NAME pyazr_ui_navigation COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ui/navigation_test.py)
set_tests_properties(pyazr_ui_navigation PROPERTIES LABELS browser SKIP_RETURN_CODE 77
                     ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# A fit continued from its written values computes the fit it continues (docs/upgrade/parameter_kinds.md section 4):
# tests/nuisance_prior_centre/continued.sh, whose part (5) needs the module (save_fit, the layouts the sessions report).
add_test(NAME pyazr_continued_fit
         COMMAND bash ${CMAKE_SOURCE_DIR}/tests/nuisance_prior_centre/continued.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(pyazr_continued_fit PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;PYAZR_PYTHON=${Python3_EXECUTABLE}"
                     PROCESSORS 2)

# The web UI's physics routes (GET /api/physics, POST /api/physics/{fn}, POST /api/session/channel_functions): the
# numbers of a direct pyazr call bit for bit, the refusals with their three parts, the bound per request.  Registered
# here, above the block that must stay last, so that it gets PYAZR_ENGINE_DIR and the skip code.
add_test(NAME pyazr_ui_physics_routes
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ui/physics_routes_test.py)
set_tests_properties(pyazr_ui_physics_routes PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
# The Physics tab's requests: every row plots without typing, the copied Python text gives the route's arrays.
add_test(NAME pyazr_ui_physics_tab
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ui/physics_tab_test.py)
set_tests_properties(pyazr_ui_physics_tab PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# The options of the calculation as <config> keys (docs/upgrade/project_format3.md section 3.12): key = keyword for every
# option and both values on tests/upgrade/T1, keyword > key > default with their sources, the refusals, set_config's
# round trip, what Project.new and save_fit write, the one table.
add_test(NAME pyazr_config_options
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/config_options_test.py)
set_tests_properties(pyazr_config_options PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The command line python -m pyazr (pyazr/cli.py): every command of the table with --json on T1/T2, describe, validate
# and schema against the program, wrong calls (exit 2 or 3 with what, allowed, fix), --help, the reference page.
add_test(NAME pyazr_cli COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cli_test.py)
set_tests_properties(pyazr_cli PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# pyazr's Objective (objective()) against today's objective().value_and_grad, residuals and objective().residuals_and_jacobian bit for
# bit on T1-T3, 13N_fast, three nuisance kinds and the gamma_decay record with priors; its rows' squares, gradient and
# Jacobian; log_prob at the sampler's excluded points and against the engine's sampler; two threads (about 55 s at one
# thread).
add_test(NAME pyazr_objective_object COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/objective_object_test.py)
set_tests_properties(pyazr_objective_object PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# The session's lock: every public member that reaches the engine takes it; two threads on one session get what each
# gets alone, bit for bit.
add_test(NAME pyazr_session_lock COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/session_lock_test.py)
set_tests_properties(pyazr_session_lock PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# The x87 precision of a thread does not reach the results: T1 in a threading.Thread whose precision is lowered to 53
# bits (as a Windows thread starts) and in the main thread lowered, bit for bit the main thread's; the control words
# of the threads printed (support/include/azure2/support/FloatingPointPrecision.h).  Two threads: the regions' scopes.
add_test(NAME pyazr_x87_precision COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/x87_precision_test.py)
set_tests_properties(pyazr_x87_precision PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77)

# pyazr's Evaluation (Session.evaluate(x)) against the program's chiSquared.out and AZUREOut_* on T1-T3, one evaluation
# per Evaluation, and the self-description of the objects.
add_test(NAME pyazr_evaluation COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/evaluation_test.py)
set_tests_properties(pyazr_evaluation PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# The three fit recipes of pyazr/examples (scipy, iminuit, emcee) on pyazr's Objective: T2 to the program's own MIGRAD
# minimum, a short emcee chain inside the support, the same with a pool of two worker processes (one thread each).
add_test(NAME pyazr_fit_recipes COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fit_recipes_test.py)
set_tests_properties(pyazr_fit_recipes PROPERTIES LABELS slow PROCESSORS 2
                     ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")

# Every script of pyazr/examples runs by its docstring's command on a project of tests/ and prints its pinned number,
# and nothing is written beside the projects (tests/pyazr/examples_test.py, one row per run): the quick rows here,
# the rows over about 4 s at one thread in pyazr_examples_slow.
add_test(NAME pyazr_examples COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/examples_test.py)
set_tests_properties(pyazr_examples PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
add_test(NAME pyazr_examples_slow COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/examples_test.py --slow)
set_tests_properties(pyazr_examples_slow PROPERTIES LABELS slow ENVIRONMENT "OMP_NUM_THREADS=1")

# The Coulomb table's build and fetch (pyazr/coulomb_table.py): a small table built, fetched by file:// and from a
# local HTTP server with its SHA-256, a corrupted copy refused and removed.  Nothing leaves this machine.
add_test(NAME pyazr_coulomb_table COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/coulomb_table_test.py)
set_tests_properties(pyazr_coulomb_table PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# pyazr as a wheel (pyproject.toml): built with pip from the source (left clean), its files against the expected
# list, unpacked (not installed) and imported alone: T2 and the stopping bit for bit as the in-tree build, the wheel
# workflow's test command (wheel_smoke.py), the wheel's AZURE2 on T1 as this build's.  Compiles the engine once more
# (two jobs: several minutes).
add_test(NAME pyazr_wheel COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wheel_test.py)
set_tests_properties(pyazr_wheel PROPERTIES LABELS slow PROCESSORS 2 TIMEOUT 7200
                     ENVIRONMENT "OMP_NUM_THREADS=1;CMAKE_BUILD_PARALLEL_LEVEL=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
# The engine's fitters on a session (docs/upgrade/fit.md): fit_with_engine against the program, fitter by fitter
# (slow), and its controls -- the fitter table, progress records, the stop request.
add_test(NAME pyazr_fit_with_engine COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fit_with_engine_test.py)
set_tests_properties(pyazr_fit_with_engine PROPERTIES LABELS slow
                     ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
# A fitted areal density in the uncertainty band (tests/layer_yield/models/fit.azr): covariance.dat names it, the band
# equals a direct propagation of the curve's differences with covariance.dat, and fixing the density narrows it.
add_test(NAME pyazr_band_target_parameters
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/band_target_parameters_test.py)
set_tests_properties(pyazr_band_target_parameters PROPERTIES SKIP_RETURN_CODE 77
                     ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
# The engine's ensemble sampler on a session (docs/upgrade/pyazr3.md section 12): sample_with_engine against the
# program's command mcmc (the chain number for number at 1 and 2 sampler threads, both bases, Park, a parameter file
# and x0, a prior of <mcmc>, a continued chain) and its controls; about 80 s.
add_test(NAME pyazr_sample_with_engine
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/sample_with_engine_test.py)
set_tests_properties(pyazr_sample_with_engine PROPERTIES PROCESSORS 2
                     ENVIRONMENT "OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>")
add_test(NAME pyazr_fit_engine_control
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fit_engine_control_test.py)
set_tests_properties(pyazr_fit_engine_control PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# A line window (capture over a broad final state's line shape): every analytic derivative of a session (gradient,
# Jacobian, model sensitivities, the fitters that read them) raises EngineError, the values keep working.
add_test(NAME pyazr_line_window_derivatives
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/line_window_derivatives_test.py)
set_tests_properties(pyazr_line_window_derivatives PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A test that cannot run its engine part -- numpy, scipy or the compiled
# _azure2 module missing (pybind11 absent at configure time builds no module),
# or no AZURE2 binary -- exits 77 after whatever pure-Python checks it has, so
# ctest lists it as Skipped rather than Passed.  An exit of 0 there once hid
# that CI had never run a single CLI == pyazr comparison.  Set on every test
# in this directory, so a new one only has to follow the convention.
# Keep this block LAST: it applies
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

# The label "long" (owner's answer of 2026-11-18, docs/upgrade/OWNER_ANSWERS.md): the tests that take the time on the
# CI's Windows and macOS runners, measured there in seconds (macOS / Windows, last CI run before 2026-11-18): thm_adjoint 2177 /
# 996, thickness 1175 / 702, session_isolation 483 / 269, ui_curves 459 / -, doppler 433 / 260, ui_utilities 278 / -,
# thm_model_average 221 / 103.  Linux runs them; Windows and macOS run `ctest -LE "slow|browser|long"`
# (.github/workflows/build.yml).  Each area keeps a test there: the THM kernels by the other pyazr_thm_* tests,
# sessions by pyazr_session_records and pyazr_session_functions, the web page by the other pyazr_ui_* tests;
# Doppler and the sample thickness by the physics regression scripts tests/doppler_fgm and tests/transmission_yield,
# which every platform runs.  A test goes in only when it takes minutes on a CI machine and its area has another test.
set_tests_properties(pyazr_thm_adjoint pyazr_thickness pyazr_session_isolation pyazr_doppler pyazr_thm_model_average
                     pyazr_ui_curves pyazr_ui_utilities PROPERTIES LABELS long)

# The words and units of pyazr's answers (the defects the example scripts met): the level scheme's pairs by key and
# kind, the units of the external-capture integrals, of an extrapolation's model values and of a gamma_decay
# parameter.
add_test(NAME pyazr_library_answers
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/library_answers_test.py)
set_tests_properties(pyazr_library_answers PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# pyazr's library calls: nothing on standard output (the engine's lines go to the logger pyazr.engine), the rate's
# files in a temporary directory of its own, a refused set-up and a wrong fit vector with their three parts.
add_test(NAME pyazr_library_calls
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/library_calls_test.py)
set_tests_properties(pyazr_library_calls PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A Project edit that leaves a phase-shift data set or extrapolation without its partial wave is refused with its
# three parts (the engine's set-up check, made at the edit), in agreement with the engine on T3.
add_test(NAME pyazr_project_phase_shift
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/project_phase_shift_test.py)
set_tests_properties(pyazr_project_phase_shift PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# A prior or a limit on a level parameter as a Project edit (<parameterSettings>): the session's objective holds the
# prior, the parameter's bounds the limits; the rows read back; replaced, removed; refusals with their three parts.
add_test(NAME pyazr_parameter_settings
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/parameter_settings_test.py)
set_tests_properties(pyazr_parameter_settings PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Settings of a potential record fitted from outside the engine (pyazr.record_fit; potential_background.md section
# 14.4): the interior charge of a wave recovered from the model's own cross sections, its interval and linear error,
# two numbers of one wave line.  Needs scipy.
add_test(NAME pyazr_record_fit
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/record_fit_test.py)
# 466 to 538 s on the macOS runner (2026-11-24): the limit, a guard against a hang, is 1800 s.
set_tests_properties(pyazr_record_fit PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" SKIP_RETURN_CODE 77 TIMEOUT 1800)
