# Taken from symengine and modified by T.P.

include(LibFindMacros)

libfind_library(gmp gmp)
set(GMP_LIBRARIES ${GMP_LIBRARY})
set(GMP_TARGETS gmp)

libfind_include(gmp.h gmp)
set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_LIBRARIES
    GMP_INCLUDE_DIRS)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
