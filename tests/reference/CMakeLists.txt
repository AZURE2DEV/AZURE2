# Reference calculation for the beam-profile experimental effect.
#
# Unlike tests/beam_profile_kernel, which pins a chi-squared produced by this
# same code, this one checks the convolution against closed-form values and
# against a second implementation of the documented integral.  It links only
# AZURE2Core and needs no model, no data files and no GUI, so it runs in
# milliseconds wherever the engine builds.

add_executable(beam_profile_reference_test beam_profile_reference_test.cpp)
target_link_libraries(beam_profile_reference_test PRIVATE AZURE2Core)

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

# CoulFuncCache: a memo key given up on keeps exact energies only, and is
# dropped when that fills without hits.
add_executable(coul_func_cache_test coul_func_cache_test.cpp)
target_link_libraries(coul_func_cache_test PRIVATE AZURE2Core)

add_test(NAME coul_func_cache COMMAND coul_func_cache_test)
set_tests_properties(coul_func_cache PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# Reference calculation for the THM distortion factor (ThmDistortion): the
# zero-range DWBA amplitude against mpmath / scipy values from
# thm_distortion_reference.py, its plane-wave limit against the closed form.
add_executable(thm_distortion_test thm_distortion_test.cpp)
target_link_libraries(thm_distortion_test PRIVATE AZURE2Core)

add_test(NAME thm_distortion COMMAND thm_distortion_test)
set_tests_properties(thm_distortion PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")

# Global optical potentials (ThmOptical): the parameters against the RIPL-3
# library and the FRONT21 code (thm_optical_reference.py), d + 40Ca elastic
# scattering against DOP853 and the measured ratio to Rutherford, and their use
# in the distortion factor.
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
# fallback of the gradient (forced by a convolved analyzing-power segment) and
# the analytic gradient agree with differences of the objective beyond the
# wall.  A copy of tests/identical_pp_res's data, adapted by the test itself.
add_executable(park_chi2_fallback_test park_chi2_fallback_test.cpp)
target_link_libraries(park_chi2_fallback_test PRIVATE AZURE2Core)
set(_park_chi2_dir ${CMAKE_CURRENT_BINARY_DIR}/park_chi2_fallback)
file(MAKE_DIRECTORY ${_park_chi2_dir}/output ${_park_chi2_dir}/checks)
file(COPY ${CMAKE_SOURCE_DIR}/tests/identical_pp_res/data DESTINATION ${_park_chi2_dir})
add_test(NAME park_chi2_fallback
         COMMAND park_chi2_fallback_test ${CMAKE_SOURCE_DIR}/tests/identical_pp_res/identical_pp_res.azr
         WORKING_DIRECTORY ${_park_chi2_dir})
set_tests_properties(park_chi2_fallback PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)

# User-defined numerical resolution function: kernel moments against closed
# forms (unit integral, centring, burst/channel widths in quadrature, energy
# interpolation) and the time-of-flight conversions.  Writes and removes a
# temporary file in its working directory.
add_executable(udr_reference_test udr_reference_test.cpp)
target_link_libraries(udr_reference_test PRIVATE AZURE2Core)
add_test(NAME udr_reference COMMAND udr_reference_test)

# Equation::Evaluate (compiled postfix program) against Equation::EvaluateTokens,
# the original string-rewriting evaluator, for every supported construct.
add_executable(equation_reference_test equation_reference_test.cpp)
target_link_libraries(equation_reference_test PRIVATE AZURE2Core)

add_test(NAME equation_reference COMMAND equation_reference_test)

# Compound stopping cross section (Bragg's rule) against the single-element
# values: per active atom and per average atom.  Reads erya/data/SRIM2013.xml,
# so it runs from external/.
if(USE_ERYA)
  add_executable(compound_stopping_reference_test compound_stopping_reference_test.cpp)
  target_include_directories(compound_stopping_reference_test PRIVATE
      ${CMAKE_SOURCE_DIR}/external/erya/include)
  target_link_libraries(compound_stopping_reference_test PRIVATE srim_utilities)
  add_test(NAME compound_stopping_reference
           COMMAND compound_stopping_reference_test
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external)
endif()

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

# CoulombTable::BuildFile (AZURE2 --build-coulomb-table, the GUI's Build Coulomb Table...):
# a small table built, counted, cancelled, and picked up by Shared() after ForgetMissing().
add_executable(coulomb_table_build_test coulomb_table_build_test.cpp)
target_link_libraries(coulomb_table_build_test PRIVATE AZURE2Core)
add_test(NAME coulomb_table_build COMMAND coulomb_table_build_test)
set_tests_properties(coulomb_table_build PROPERTIES TIMEOUT 120)

# AZUREParams::ReadUserParameters on a param.sav with "nan" errors (written after a
# fit with an invalid covariance), blank and malformed lines: it used to loop
# until std::bad_alloc.  The timeout turns a hang into a failure.
add_executable(parameter_file_reference_test parameter_file_reference_test.cpp)
target_link_libraries(parameter_file_reference_test PRIVATE AZURE2Core)
add_test(NAME parameter_file_reference COMMAND parameter_file_reference_test)
set_tests_properties(parameter_file_reference PROPERTIES TIMEOUT 10)
