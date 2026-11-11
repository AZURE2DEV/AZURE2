# stopping/: AZURE2's stopping and straggling (docs/upgrade/stopping.md, docs/upgrade/ARCHITECTURE.md section 7).
# A leaf library: it uses physics/ (constants, element symbols) and the C++ standard library, nothing else of AZURE2
# (repo_consistency check 32 refuses any other include).  The engine (AZURE2Core) links it.
#
# Layout: public headers in include/azure2/stopping/, sources in src/, SRIM-2013's data files in data/srim/ (compiled
# in at build time, so the library needs no file at run time), tests in tests/, the command-line tool and the
# generator of the LunaSRIM reference values in tools/.

# SRIM's data files, compiled in: TableText.cpp is generated from data/srim/ (rebuilt when a file changes).  Each file
# becomes a list of raw string literals of at most 30000 characters (compilers limit one literal), carriage returns
# removed.  Without the files (the owner may decide not to ship them) the build embeds nothing and the library reads
# them from the directory in the environment variable AZURE2_SRIM_DATA_DIR (stopping/data/srim/README.md).
function(azure2_stopping_pieces file out_var)
  set(pieces "")
  if(EXISTS ${file})
    file(READ ${file} text)
    string(REPLACE "\r" "" text "${text}")
    string(LENGTH "${text}" length)
    set(pos 0)
    while(pos LESS length)
      string(SUBSTRING "${text}" ${pos} 30000 piece)
      string(APPEND pieces "    R\"azure2srim(${piece})azure2srim\",\n")
      math(EXPR pos "${pos} + 30000")
    endwhile()
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${file})
  endif()
  set(${out_var} "${pieces}" PARENT_SCOPE)
endfunction()
set(AZURE2_STOPPING_SRIM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data/srim)
azure2_stopping_pieces(${AZURE2_STOPPING_SRIM_DIR}/SCOEF03.DAT AZURE2_STOPPING_SCOEF_PIECES)
azure2_stopping_pieces(${AZURE2_STOPPING_SRIM_DIR}/SNUC03.dat AZURE2_STOPPING_SNUC_PIECES)
if(AZURE2_STOPPING_SCOEF_PIECES STREQUAL "" OR AZURE2_STOPPING_SNUC_PIECES STREQUAL "")
  message(STATUS "stopping: SRIM-2013 data files not found in ${AZURE2_STOPPING_SRIM_DIR}; the library will read them "
                 "from AZURE2_SRIM_DATA_DIR")
  set(AZURE2_STOPPING_SCOEF_PIECES "")
  set(AZURE2_STOPPING_SNUC_PIECES "")
endif()
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

# The tests (sources named in full: repo_consistency checks that every tests/*_test.cpp is named here).  They read
# tests/data/srim_reference.csv (tools/make_srim_reference.py: LunaSRIM's own values) and need no Python.
#   stopping_tables          SRIM's data files: parsed completely, read from a directory as from the compiled copy
#   stopping_srim_directory  the same test with AZURE2_SRIM_DATA_DIR set: SrimTables() reads the directory
#   stopping_electronic      electronic stopping against LunaSRIM, every branch and boundary, phase factor
#   stopping_nuclear         SNUC03.dat nuclear stopping against LunaSRIM, its definition, the form change at 30
#   stopping_material        formulas, default phase, compounds against LunaSRIM, uncertainty, user tables, provenance
#   stopping_energy_loss     N, E_out, range, scale derivatives (closed forms), break points, layers against LunaSRIM
#   stopping_straggling      Bohr, TRIM's model, Tschalaer against closed forms, layers against LunaSRIM
set(AZURE2_STOPPING_REFERENCE ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/srim_reference.csv)
foreach(_src tables_test.cpp nuclear_test.cpp electronic_test.cpp material_test.cpp energy_loss_test.cpp
             straggling_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(stopping_${_t}_test tests/${_src})
  target_link_libraries(stopping_${_t}_test PRIVATE azure2_stopping)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(stopping_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME stopping_${_t} COMMAND stopping_${_t}_test ${AZURE2_STOPPING_REFERENCE} ${AZURE2_STOPPING_SRIM_DIR})
endforeach()
add_test(NAME stopping_srim_directory COMMAND stopping_tables_test ${AZURE2_STOPPING_REFERENCE} ${AZURE2_STOPPING_SRIM_DIR})
set_tests_properties(stopping_srim_directory PROPERTIES ENVIRONMENT "AZURE2_SRIM_DATA_DIR=${AZURE2_STOPPING_SRIM_DIR}")

# The command-line tool: epsilon, energy loss, range and straggling tables (docs/upgrade/stopping.md section 11).
add_executable(azure2-stopping tools/azure2_stopping.cpp)
target_link_libraries(azure2-stopping PRIVATE azure2_stopping)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2-stopping PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME stopping_cli COMMAND azure2-stopping stopping --projectile a --material Ta2O5 --from 0.1 --to 1 --points 3)
