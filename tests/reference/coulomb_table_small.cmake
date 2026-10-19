# ctest driver for coulomb_table_small: build a small table, then validate it.
# 30 nodes per decade gives ~1e-6 interpolation errors over this range; the
# limit only has to catch a broken build, lookup or file round trip.
execute_process(COMMAND ${TOOL} build ${OUT} --lmax 2 --eta 0.1 3 --rho 0.1 10 --per-decade 30 --threads 2
                RESULT_VARIABLE rb)
if(NOT rb EQUAL 0)
  message(FATAL_ERROR "azure2-coulomb-table build failed (${rb})")
endif()
execute_process(COMMAND ${TOOL} validate ${OUT} --n 4000 --limit 1e-4 --threads 2 RESULT_VARIABLE rv)
file(REMOVE ${OUT})
if(NOT rv EQUAL 0)
  message(FATAL_ERROR "azure2-coulomb-table validate failed (${rv})")
endif()
