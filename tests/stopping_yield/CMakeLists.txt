# tests/stopping_yield (README.md): the convention of the stopping per active nucleus, against numbers put together by
# hand from the library's per-atom values.  check.sh (tests/run_tests.sh) runs the projects.
add_executable(stopping_convention_test stopping_convention_test.cpp)
target_link_libraries(stopping_convention_test PRIVATE AZURE2Core)
add_test(NAME stopping_convention COMMAND stopping_convention_test)
