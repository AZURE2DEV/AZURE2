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
