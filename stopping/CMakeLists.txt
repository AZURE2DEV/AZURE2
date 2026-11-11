# stopping/: AZURE2's stopping and straggling (docs/upgrade/stopping.md, docs/upgrade/ARCHITECTURE.md section 7).
# A leaf library: it uses physics/ (constants, element symbols) and the C++ standard library, nothing else of AZURE2
# (repo_consistency check 32 refuses any other include).  Nothing in the engine links it yet.
#
# Layout: public headers in include/azure2/stopping/, sources in src/, the data tables in data/ (compiled in at build
# time from the .tsv files, so the tables are the files), tests in tests/, the command-line tool and the comparison
# with ERYA in tools/.

# The tables, compiled in: TableText.cpp is generated from data/*.tsv (rebuilt when a table changes).
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/data/elements.tsv AZURE2_STOPPING_ELEMENTS_TSV)
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/data/electronic_coefficients.tsv AZURE2_STOPPING_COEFFICIENTS_TSV)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
             ${CMAKE_CURRENT_SOURCE_DIR}/data/elements.tsv ${CMAKE_CURRENT_SOURCE_DIR}/data/electronic_coefficients.tsv)
configure_file(src/TableText.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/generated/TableText.cpp @ONLY)

add_library(azure2_stopping STATIC
  src/Electronic.cpp
  src/EnergyLoss.cpp
  src/Material.cpp
  src/Nuclear.cpp
  src/Stopping.cpp
  src/Straggling.cpp
  src/Tables.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/generated/TableText.cpp)
add_library(azure2::stopping ALIAS azure2_stopping)
target_include_directories(azure2_stopping PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(azure2_stopping PUBLIC azure2_physics)
target_compile_features(azure2_stopping PUBLIC cxx_std_17)
set_target_properties(azure2_stopping PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_stopping PRIVATE -Wall -Wextra -Wpedantic)
endif()

# The tests (sources named in full: repo_consistency checks that every tests/*_test.cpp is named here).
#   stopping_tables       the compiled tables: 92 elements, 184 rows, sources, parser refusals
#   stopping_nuclear      ZBL fit against the universal potential's scattering integral and an independent evaluation
#   stopping_electronic   fit forms, Bethe terms (Bloch against the digamma function), joins, velocity scaling
#   stopping_material     formulas, Bragg additivity, uncertainty, user tables and formulas, provenance
#   stopping_energy_loss  N, E_out, range and scale derivatives against closed forms for epsilon = c E^p
#   stopping_straggling   Bohr, the thin-layer and constant-S limits, Tschalaer against a closed form, Chu, Yang
foreach(_src tables_test.cpp nuclear_test.cpp electronic_test.cpp material_test.cpp energy_loss_test.cpp
             straggling_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(stopping_${_t}_test tests/${_src})
  target_link_libraries(stopping_${_t}_test PRIVATE azure2_stopping)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(stopping_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME stopping_${_t} COMMAND stopping_${_t}_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/reference.csv)
endforeach()

# The command-line tool: epsilon, energy loss, range and straggling tables (docs/upgrade/stopping.md section 11).
add_executable(azure2-stopping tools/azure2_stopping.cpp)
target_link_libraries(azure2-stopping PRIVATE azure2_stopping)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2-stopping PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME stopping_cli COMMAND azure2-stopping stopping --projectile a --material Ta2O5 --from 0.1 --to 1 --points 3)

# The comparison with ERYA (stopping.md section 9): built only where ERYA is (USE_ERYA), run by hand with
# tools/compare_erya.sh; not a test.
if(TARGET srim_utilities)
  add_executable(stopping_erya_compare tools/erya_compare.cpp)
  target_link_libraries(stopping_erya_compare PRIVATE azure2_stopping srim_utilities)
  target_include_directories(stopping_erya_compare PRIVATE ${CMAKE_SOURCE_DIR}/external/erya/include)
endif()
