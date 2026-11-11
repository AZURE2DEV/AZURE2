add_library(coul complex_functions.cpp cwfcomp.cpp cwfcomp_accurate.cpp)
# cwfcomp_accurate.cpp calls GSL's Bessel functions (heavy-ion threshold fallback, dev 56512d6).  Declaring the
# dependency here makes CMake put GSL AFTER libcoul.a on every link line; without it GNU ld on a static GSL (the
# MinGW/Windows CI) left gsl_sf_bessel_I0/I1/K0/K1 unresolved, because AZURE2Core lists GSL before coul.
target_link_libraries(coul PUBLIC ${GSL_LIBRARIES})
