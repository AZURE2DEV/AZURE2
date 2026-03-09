#It looks like a GSL environment variable is now available GSL_DIR, so I am using this to elimiate the hard coded path

find_path(GSL_INCLUDE_DIR gsl/gsl_cblas.h
          HINTS $ENV{HOME}/local/include $ENV{HOME}/include $ENV{GSL_DIR}/include) # /opt/crc/scilib/gsl/1.15/gnu-4.6.2/include
# /opt/crc/scilib/gsl/1.15/gnu-4.6.2/include

find_library(GSL_LIBRARY NAMES gsl 
             HINTS $ENV{HOME}/local/lib $ENV{HOME}/lib $ENV{GSL_DIR}/lib) # /opt/crc/scilib/gsl/1.15/gnu-4.6.2/lib
# /opt/crc/scilib/gsl/1.15/gnu-4.6.2/lib

find_library(GSLCBLAS_LIBRARY NAMES gslcblas 
             HINTS $ENV{HOME}/local/lib $ENV{HOME}/lib $ENV{GSL_DIR}/lib) # /opt/crc/scilib/gsl/1.15/gnu-4.6.2/lib
# /opt/crc/scilib/gsl/1.15/gnu-4.6.2/lib

set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(gsl  DEFAULT_MSG
                                  GSL_LIBRARY GSL_INCLUDE_DIR)

mark_as_advanced(GSL_INCLUDE_DIR GSL_LIBRARY GSLCBLAS_LIBRARY)
