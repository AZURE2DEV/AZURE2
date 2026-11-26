# Reference calculation for the beam-profile experimental effect.
#
# Unlike tests/beam_profile_kernel, which pins a chi-squared produced by this
# same code, this one checks the convolution against closed-form values and
# against a second implementation of the documented integral.  It links only
# AZURE2Core and needs no model, no data files and no GUI, so it runs in
# milliseconds wherever the engine builds.

add_executable(beam_profile_reference_test beam_profile_reference_test.cpp)
target_link_libraries(beam_profile_reference_test PRIVATE AZURE2Core azure2_convert_specs)

add_test(NAME beam_profile_reference COMMAND beam_profile_reference_test)

# Reference calculation for the THM external Coulomb term (ThmCoulombTerm),
# against mpmath values from thm_coulomb_term_reference.py.
add_executable(thm_coulomb_term_test thm_coulomb_term_test.cpp)
target_link_libraries(thm_coulomb_term_test PRIVATE AZURE2Core)

add_test(NAME thm_coulomb_term COMMAND thm_coulomb_term_test)

# The shared threshold helper ChannelFunc: S, P and dS/dE of a channel at,
# just below and just above threshold (finite, continuous, fast), and bit for
# bit the old routines away from it.
add_executable(channel_threshold_test channel_threshold_test.cpp)
target_link_libraries(channel_threshold_test PRIVATE AZURE2Core)

add_test(NAME channel_threshold COMMAND channel_threshold_test)

# GslQuiet: quiet GSL calls from many OpenMP threads leave the process-wide GSL
# error handler as it was (a bare save/off/restore could lose it for good).
add_executable(gsl_quiet_test gsl_quiet_test.cpp)
target_link_libraries(gsl_quiet_test PRIVATE AZURE2Core)

add_test(NAME gsl_quiet COMMAND gsl_quiet_test)
set_tests_properties(gsl_quiet PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")

# Reference calculation for the THM distortion factor (ThmDistortion): the
# zero-range DWBA amplitude against mpmath / scipy values from
# thm_distortion_reference.py, its plane-wave limit against the closed form.
add_executable(thm_distortion_test thm_distortion_test.cpp)
target_link_libraries(thm_distortion_test PRIVATE AZURE2Core)

add_test(NAME thm_distortion COMMAND thm_distortion_test)
set_tests_properties(thm_distortion PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Global optical potentials (physics/transfer/GlobalOptical.h; the parameters
# against RIPL-3 and FRONT21 are physics_transfer_global_optical): d + 40Ca
# elastic scattering against DOP853 and the measured ratio to Rutherford, and
# their use in the distortion factor.
# The refusals of the THM set-up name the way on (what is wrong, what is allowed, how to go on).
add_executable(thm_refusals_test thm_refusals_test.cpp)
target_link_libraries(thm_refusals_test PRIVATE AZURE2Core)
add_test(NAME thm_refusals COMMAND thm_refusals_test)

add_executable(thm_optical_test thm_optical_test.cpp)
target_link_libraries(thm_optical_test PRIVATE AZURE2Core)

add_test(NAME thm_optical COMMAND thm_optical_test)
set_tests_properties(thm_optical PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Reference calculation for the fixed-angle THM observable (ThmAngular): the
# Legendre form and its window average against the literal M-sum of
# thm_fixed_angle_reference.py (sympy Clebsch-Gordan, mpmath quadrature).
add_executable(thm_fixed_angle_test thm_fixed_angle_test.cpp)
target_link_libraries(thm_fixed_angle_test PRIVATE AZURE2Core)

add_test(NAME thm_fixed_angle COMMAND thm_fixed_angle_test)

# Reference calculation for the distorted-wave THM entrance vertex
# (ThmDwVertex): its plane-wave limit against the plane-wave vertex, and the
# Coulomb-distorted Gram matrix against a direct quadrature of the source
# (thm_dw_vertex_reference.py).
add_executable(thm_dw_vertex_test thm_dw_vertex_test.cpp)
target_link_libraries(thm_dw_vertex_test PRIVATE AZURE2Core)

add_test(NAME thm_dw_vertex COMMAND thm_dw_vertex_test)
set_tests_properties(thm_dw_vertex PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# The spectator-direction window of a THM experiment (spectatorAngles=,
# ThmDistortion::AngleNodes): the measure of the accepted c.m. directions
# against a brute-force sum, lab windows with one and two c.m. branches, the q
# cut of a ps window, zero-width windows and acceptance tables.
add_executable(thm_spectator_angles_test thm_spectator_angles_test.cpp)
target_link_libraries(thm_spectator_angles_test PRIVATE AZURE2Core)

add_test(NAME thm_spectator_angles COMMAND thm_spectator_angles_test)

# The Park J > 0 penalty in AZURECalc::Chi2Value: the finite-difference
# fallback of the gradient (taken for a folded capture analyzing power,
# tests/13N_capture_ay with a beam spread) and the analytic gradient (p+p,
# plain and with a convolved analyzing power) agree with differences of the
# objective beyond the wall; and the first project without Park under --use-rmc, whose adjoint bails at a gate of
# AMatrixFunc::PointTBar.  Copies of tests/identical_pp_res's and tests/13N_capture_ay's data, adapted by the test
# itself.
add_executable(park_chi2_fallback_test park_chi2_fallback_test.cpp)
target_link_libraries(park_chi2_fallback_test PRIVATE AZURE2Core)
set(_park_chi2_dir ${CMAKE_CURRENT_BINARY_DIR}/park_chi2_fallback)
file(MAKE_DIRECTORY ${_park_chi2_dir}/output ${_park_chi2_dir}/checks)
file(COPY ${CMAKE_SOURCE_DIR}/tests/identical_pp_res/data DESTINATION ${_park_chi2_dir})
file(COPY ${CMAKE_SOURCE_DIR}/tests/13N_capture_ay/data DESTINATION ${_park_chi2_dir})
add_test(NAME park_chi2_fallback
         COMMAND park_chi2_fallback_test ${CMAKE_SOURCE_DIR}/tests/identical_pp_res/identical_pp_res.azr
                 ${CMAKE_SOURCE_DIR}/tests/13N_capture_ay/13N_capture_ay.azr
         WORKING_DIRECTORY ${_park_chi2_dir})
set_tests_properties(park_chi2_fallback PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)

# A non-finite parameter vector (a fit step to NaN or infinity) is a failed evaluation, never an abort: the
# chi-squared ceiling, NaN, false; a parameter file with a non-finite value is refused.  A synthetic model with a
# potential background (tests/potential_background/models/fit_level.azr) and the same without the potential.
add_executable(nonfinite_step_test nonfinite_step_test.cpp)
target_link_libraries(nonfinite_step_test PRIVATE AZURE2Core)
set(_nonfinite_dir ${CMAKE_CURRENT_BINARY_DIR}/nonfinite_step)
file(MAKE_DIRECTORY ${_nonfinite_dir}/output ${_nonfinite_dir}/checks)
file(COPY ${CMAKE_SOURCE_DIR}/tests/potential_background/models/data DESTINATION ${_nonfinite_dir})
add_test(NAME nonfinite_step
         COMMAND nonfinite_step_test ${CMAKE_SOURCE_DIR}/tests/potential_background/models/fit_level.azr
         WORKING_DIRECTORY ${_nonfinite_dir})
set_tests_properties(nonfinite_step PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)

# A shift column's difference step stays inside the planned nodes (ShiftedFold.h, question 5): a varied shift without
# a prior on a thin layer (tests/thin_layer_limit's model at three densities), its gradient and residual Jacobian near
# the edge of the reach against central differences.  Before, the process ended (std::logic_error in the refold).
add_executable(shift_column_inside_test shift_column_inside_test.cpp)
target_link_libraries(shift_column_inside_test PRIVATE AZURE2Core)
set(_shift_column_dir ${CMAKE_CURRENT_BINARY_DIR}/shift_column_inside)
file(MAKE_DIRECTORY ${_shift_column_dir}/output ${_shift_column_dir}/checks)
file(COPY ${CMAKE_SOURCE_DIR}/tests/thin_layer_limit/models/data DESTINATION ${_shift_column_dir})
add_test(NAME shift_column_inside
         COMMAND shift_column_inside_test ${CMAKE_SOURCE_DIR}/tests/thin_layer_limit/models/thin_limit.azr
         WORKING_DIRECTORY ${_shift_column_dir})
set_tests_properties(shift_column_inside PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)

# User-defined numerical resolution function: kernel moments against closed
# forms (unit integral, centring, burst/channel widths in quadrature, energy
# interpolation) and the time-of-flight conversions.  Writes and removes a
# temporary file in its working directory.
# A <targetInt> line that ends inside its first nine fields: every TargetEffect member has its initialised value
# (step 3a; before, they were left indeterminate), and a complete line still sets them.
add_executable(target_effect_line_test target_effect_line_test.cpp)
target_link_libraries(target_effect_line_test PRIVATE AZURE2Core azure2_convert_specs)
add_test(NAME target_effect_line COMMAND target_effect_line_test)

# The old-line adaptor (convert/'s SpecsFromOldLines.h) and Segment's spec constructors: the isDiff table, components,
# the reading refusals, prior centres, abundances, <targetInt> data-set lists.
add_executable(old_line_specs_test old_line_specs_test.cpp)
target_link_libraries(old_line_specs_test PRIVATE AZURE2Core azure2_convert_specs)
add_test(NAME old_line_specs COMMAND old_line_specs_test)

add_executable(udr_reference_test udr_reference_test.cpp)
target_link_libraries(udr_reference_test PRIVATE AZURE2Core azure2_convert_specs)
add_test(NAME udr_reference COMMAND udr_reference_test)

# Equation::Evaluate (compiled postfix program) against Equation::EvaluateTokens,
# the original string-rewriting evaluator, for every supported construct.
add_executable(equation_reference_test equation_reference_test.cpp)
target_link_libraries(equation_reference_test PRIVATE AZURE2Core)

add_test(NAME equation_reference COMMAND equation_reference_test)

# Closed-channel shift function of neutral (eta = 0) pairs,
# ShftFunc::NeutralBoundShift, against the numerical Whittaker derivative it
# replaces and the l = 0, 1 textbook expressions; prints the speed-up.
add_executable(neutral_channel_reference_test neutral_channel_reference_test.cpp)
target_link_libraries(neutral_channel_reference_test PRIVATE AZURE2Core)
add_test(NAME neutral_channel_reference COMMAND neutral_channel_reference_test)

# Universal Coulomb tables (CoulombTable, azure2-coulomb-table): a small table
# (l <= 2, eta in [0.1, 3], rho in [0.1, 10], 30 nodes per decade) built and validated in
# a few seconds; the full-size validation is in docs/source/developer/coulomb_table.rst.
add_test(NAME coulomb_table_small
         COMMAND ${CMAKE_COMMAND} -DTOOL=$<TARGET_FILE:azure2-coulomb-table>
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/coulomb_table_small.bin
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/coulomb_table_small.cmake)

# CoulombTable::BuildFile (AZURE2 --build-coulomb-table, the setup program's Run > Build the Coulomb Table...):
# a small table built, counted, cancelled, and picked up by Shared() after ForgetMissing().
add_executable(coulomb_table_build_test coulomb_table_build_test.cpp)
target_link_libraries(coulomb_table_build_test PRIVATE AZURE2Core)
add_test(NAME coulomb_table_build COMMAND coulomb_table_build_test)
set_tests_properties(coulomb_table_build PROPERTIES TIMEOUT 120)

# CoulombTable::LookupOpen deep below the barrier: the table's phase (stored as phi - rho) is rounding noise where
# phi < 1e-16, so below CoulombTable::kSmallPhase the lookup takes it from CoulombTable::SmallPhase (power series).
add_executable(coulomb_table_small_phase_test coulomb_table_small_phase_test.cpp)
target_link_libraries(coulomb_table_small_phase_test PRIVATE AZURE2Core)
add_test(NAME coulomb_table_small_phase COMMAND coulomb_table_small_phase_test)
set_tests_properties(coulomb_table_small_phase PROPERTIES TIMEOUT 300)

# AZUREParams::ReadUserParameters on a param.sav with "nan" errors (written after a
# fit with an invalid covariance), blank and malformed lines: it used to loop
# until std::bad_alloc.  The timeout turns a hang into a failure.
add_executable(parameter_file_reference_test parameter_file_reference_test.cpp)
target_link_libraries(parameter_file_reference_test PRIVATE AZURE2Core)
add_test(NAME parameter_file_reference COMMAND parameter_file_reference_test)
set_tests_properties(parameter_file_reference PROPERTIES TIMEOUT 10)

# The adaptive sub-point grid of the experimental effects (AdaptiveIntegrationGrid):
# a Lorentzian moved across a fixed window, the grid laid out again at each position
# as WorkingCopy::ReplaceSubPointGrids does, its trapezoid integral against the closed form.
# Checks the accuracy and that the error does not swing as the level moves.
add_executable(adaptive_grid_reference_test adaptive_grid_reference_test.cpp)
target_link_libraries(adaptive_grid_reference_test PRIVATE AZURE2Core)
add_test(NAME adaptive_grid_reference COMMAND adaptive_grid_reference_test)

# Wigner limit of a channel: the particle formula, and no limit for photon, beta and
# radius-less channels (a photon pair's mu = 0 gave inf before 2026-10-27, audit item 32).
add_executable(wigner_limit_reference_test wigner_limit_reference_test.cpp)
target_link_libraries(wigner_limit_reference_test PRIVATE AZURE2Core)
add_test(NAME wigner_limit_reference COMMAND wigner_limit_reference_test)

# The integrated Fermi function stays finite down to the beta end point (it was NaN within ~0.2 keV of it).
add_executable(fermi_endpoint_reference_test fermi_endpoint_reference_test.cpp)
target_link_libraries(fermi_endpoint_reference_test PRIVATE AZURE2Core)
add_test(NAME fermi_endpoint_reference COMMAND fermi_endpoint_reference_test)

# InputFiles' cache must not serve a same-size rewrite whose stamp did not change (the Windows CI failure of
# pyazr_natural_target, 2026-11-05), and must still serve an old, unchanged file from memory.
add_executable(input_files_cache_test input_files_cache_test.cpp)
target_link_libraries(input_files_cache_test PRIVATE AZURE2Core)
add_test(NAME input_files_cache COMMAND input_files_cache_test)


# physics/coulomb's GSL-free Whittaker function against the engine's GSL path (gsl_sf_hyperg_U_e10_e) on the mpmath
# grid of physics/tests/data/whittaker.csv: agreement within GSL's own error; prints both errors and the time per
# call (docs/upgrade/coulomb.md section 1.4).
add_executable(whittaker_gsl_test whittaker_gsl_test.cpp)
target_link_libraries(whittaker_gsl_test PRIVATE AZURE2Core)
add_test(NAME whittaker_gsl COMMAND whittaker_gsl_test ${CMAKE_SOURCE_DIR}/physics/tests/data/whittaker.csv)

# The ensemble sampler samples the fit's objective and owns its state (docs/upgrade/objective.md section 6): on the
# gamma_decay record project of tests/gamma_decay (with the record's priors), in Brune's and Park's bases and both
# sampling bases, ln P = -chi^2/2 + C at three points and -infinity outside the record's limits; a short chain stays
# inside them; the record's priors move ln P by -1/2 x their term; T2's identity with its Gaussian priors' constant;
# two samplers in two threads write the chains of two runs one after the other; the support is the fit's limits of
# every kind (ln P = -chi^2/2 + C just inside, -infinity just outside, both bases); a chain started next to a limit.
if(USE_MCMC)
  add_executable(sampler_objective_test sampler_objective_test.cpp)
  target_link_libraries(sampler_objective_test PRIVATE AZURE2Core)
  set(_sampler_objective_dir ${CMAKE_CURRENT_BINARY_DIR}/sampler_objective)
  file(MAKE_DIRECTORY ${_sampler_objective_dir}/output ${_sampler_objective_dir}/checks)
  file(COPY ${CMAKE_SOURCE_DIR}/tests/gamma_decay/models/data DESTINATION ${_sampler_objective_dir})
  file(COPY ${CMAKE_SOURCE_DIR}/tests/upgrade/T2/data DESTINATION ${_sampler_objective_dir})
  foreach(_data 13N/data screening_6li/data stopping_yield/models/data rectangular_resolution/models/data)
    file(COPY ${CMAKE_SOURCE_DIR}/tests/${_data} DESTINATION ${_sampler_objective_dir})
  endforeach()
  file(COPY ${CMAKE_SOURCE_DIR}/tests/upgrade/T2/cov_diff.txt ${CMAKE_SOURCE_DIR}/tests/upgrade/T2/prior_cov.f3.txt
       DESTINATION ${_sampler_objective_dir})
  add_test(NAME sampler_objective
           COMMAND sampler_objective_test ${CMAKE_SOURCE_DIR}/tests/gamma_decay/models/record.azr
                   ${CMAKE_SOURCE_DIR}/tests/upgrade/T2/T2.azr
                   ${CMAKE_SOURCE_DIR}/tests/mcmc_physical_norm/13N_fast.azr
                   ${CMAKE_SOURCE_DIR}/tests/screening_6li/fit.azr
                   ${CMAKE_SOURCE_DIR}/tests/stopping_yield/models/fit.azr
                   ${CMAKE_SOURCE_DIR}/tests/rectangular_resolution/models/fit.azr
           WORKING_DIRECTORY ${_sampler_objective_dir})
  # 248 to 526 s on the macOS runner (five runs of 2026-11-24): the limit, a guard against a hang, is 1800 s.
  set_tests_properties(sampler_objective PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 1800)
endif()

# MigradSeed: the bounded curvature of MIGRAD's seed against hand-computed values, and the watch of the first step
# (a failed first line search is recognised and stops the pass) as a state machine.
add_executable(migrad_seed_test migrad_seed_test.cpp)
target_link_libraries(migrad_seed_test PRIVATE AZURE2Core)

add_test(NAME migrad_seed COMMAND migrad_seed_test)
