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
