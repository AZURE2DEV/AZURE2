find_path(MINUIT2_INCLUDE_DIR Minuit2/MnUserFcn.h
          HINTS $ENV{HOME}/local/include $ENV{HOME}/include $ENV{ROOTSYS}/include/root) #/afs/crc.nd.edu/user/n/nsl/nuclear/x86_64_linux_el6/install/root_v5.34.19/include/root

find_library(MINUIT2_LIBRARY NAMES Minuit2
             HINTS $ENV{HOME}/local/lib $ENV{HOME}/lib $ENV{ROOTSYS}/lib/root) #/afs/crc.nd.edu/user/n/nsl/nuclear/x86_64_linux_el6/install/root_v5.34.19/lib/root

set (MINUIT2_LIBRARIES ${MINUIT2_LIBRARY})

if (MINUIT2_INCLUDE_DIR MATCHES $ENV{ROOTSYS}/include)
execute_process(COMMAND root-config --libs 
    RESULT_VARIABLE return_code
    OUTPUT_VARIABLE output ERROR_VARIABLE error
    OUTPUT_STRIP_TRAILING_WHITESPACE)
if (output)
   set (MINUIT2_LIBRARIES ${MINUIT2_LIBRARIES} ${output})
endif (output)
endif (MINUIT2_INCLUDE_DIR MATCHES $ENV{ROOTSYS}/include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Minuit2  DEFAULT_MSG
                                  MINUIT2_LIBRARY MINUIT2_INCLUDE_DIR)

mark_as_advanced(MINUIT2_INCLUDE_DIR MINUIT2_LIBRARY)
