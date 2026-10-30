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

# calculate_chi2 (physical parameters) is the same objective as calculate_chi2_rwa:
# the normalizations in its vector are applied (they were not before 2026-10-26),
# the Park / abundance / scattering-length terms are included, a rejected vector is +inf.
add_test(NAME pyazr_physical_norm
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/physical_norm_test.py)
set_tests_properties(pyazr_physical_norm PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

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
set_tests_properties(pyazr_inactive_lines PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

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
set_tests_properties(pyazr_output_files PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The sqrt(E) energy-shift block on a segment line: read, written, and (with
# the engine) counted and penalised exactly as AZURE2 does.
add_test(NAME pyazr_sqrt_shift
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/sqrt_shift_test.py)
set_tests_properties(pyazr_sqrt_shift PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

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
set_tests_properties(pyazr_park_gradient PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")

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
set_tests_properties(pyazr_sav_readers PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

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

# chiSquared.out written through the API carries the Park, abundance and scattering-length totals
# (2026-10-27, audit item 30).
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
