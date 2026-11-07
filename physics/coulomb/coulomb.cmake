# physics/coulomb: the Coulomb functions F, F', G, G' and sigma_l (docs/upgrade/coulomb.md), included into the
# physics target by physics/CMakeLists.txt (include(), so that the per-source properties below reach the target).
# The third-party sources sit in one directory per algorithm (coulsm/, coul/, coulfg/, coulcc/), taken from pycoulomb
# (README.md); coul/ and coulsm/ have headers of the same names, so each algorithm's sources see only its own include
# directory.  The public header is include/azure2/physics/coulomb/Coulomb.h.
set(_pc ${CMAKE_CURRENT_LIST_DIR})
set(_coulsm ${_pc}/coulsm/src/complex_functions.cpp ${_pc}/coulsm/src/cwfcomp.cpp
            ${_pc}/coulsm/src/cwfcomp_accurate.cpp ${_pc}/src/method_coulsm.cpp)
set(_coul ${_pc}/coul/src/complex_functions.cpp ${_pc}/coul/src/cwfcomp.cpp ${_pc}/coul/src/ode_int.cpp
          ${_pc}/src/method_coul.cpp)
set(_coulfg ${_pc}/coulfg/src/coulfg2.c ${_pc}/src/method_coulfg.cpp)
set(_coulcc ${_pc}/coulcc/src/coulcc.cpp ${_pc}/coulcc/src/kernels.cpp ${_pc}/src/method_coulcc.cpp)
target_sources(azure2_physics PRIVATE ${_pc}/src/Coulomb.cpp ${_coulsm} ${_coul} ${_coulfg} ${_coulcc})
target_include_directories(azure2_physics PUBLIC ${_pc}/include)
# Include directories: the third-party sources see their own directory; the adaptors (ours, compiled with the physics
# target's warnings) see it as a system directory (-isystem), so that the third-party headers' warnings are not
# reported against them.
foreach(_m coulsm coul coulfg coulcc)
  set(_dirs ${_pc}/${_m}/include)
  if(_m STREQUAL "coulcc")
    list(APPEND _dirs ${_pc}/coulcc/src)
  endif()
  set(_srcs ${_${_m}})
  list(FILTER _srcs EXCLUDE REGEX "/src/method_")
  set_source_files_properties(${_srcs} PROPERTIES INCLUDE_DIRECTORIES "${_dirs}")
  set_source_files_properties(${_pc}/src/method_${_m}.cpp PROPERTIES INCLUDE_DIRECTORIES "${_pc}/src;${_pc}/include")
  foreach(_d ${_dirs})
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      set_property(SOURCE ${_pc}/src/method_${_m}.cpp APPEND PROPERTY COMPILE_OPTIONS -isystem ${_d})
    else()
      set_property(SOURCE ${_pc}/src/method_${_m}.cpp APPEND PROPERTY INCLUDE_DIRECTORIES ${_d})
    endif()
  endforeach()
endforeach()
# The third-party code is kept as published: its warnings are not ours to fix (the adaptors and Coulomb.cpp are
# compiled with the physics target's -Wall -Wextra -Wpedantic).  COULCC reproduces the Fortran statement by
# statement only without fused multiply-add contraction (coulcc/ README of pycoulomb's branch coulcc-cpp).
set(_third_party ${_pc}/coulsm/src/complex_functions.cpp ${_pc}/coulsm/src/cwfcomp.cpp ${_pc}/coulsm/src/cwfcomp_accurate.cpp
                 ${_pc}/coul/src/complex_functions.cpp ${_pc}/coul/src/cwfcomp.cpp ${_pc}/coul/src/ode_int.cpp
                 ${_pc}/coulfg/src/coulfg2.c ${_pc}/coulcc/src/coulcc.cpp ${_pc}/coulcc/src/kernels.cpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set_property(SOURCE ${_third_party} APPEND PROPERTY COMPILE_OPTIONS -w)
  set_property(SOURCE ${_pc}/coulcc/src/coulcc.cpp ${_pc}/coulcc/src/kernels.cpp APPEND PROPERTY COMPILE_OPTIONS
               -ffp-contract=off -fno-fast-math)
endif()
if(WIN32)
  # M_PI and friends under -std=c++17 with MinGW-w64 and MSVC.
  set_property(SOURCE ${_third_party} ${_pc}/src/method_coulsm.cpp ${_pc}/src/method_coul.cpp APPEND PROPERTY
               COMPILE_DEFINITIONS _USE_MATH_DEFINES)
endif()
