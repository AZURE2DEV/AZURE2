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
