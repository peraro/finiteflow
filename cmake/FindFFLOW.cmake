# Taken from symengine

include(LibFindMacros)

libfind_include(fflow/capi_native_ext.h fflow)
libfind_library(fflow fflow)

set(FFLOW_LIBRARIES ${FFLOW_LIBRARY})
set(FFLOW_INCLUDE_DIRS ${FFLOW_INCLUDE_DIR})
set(FFLOW_TARGETS fflow)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFLOW DEFAULT_MSG FFLOW_LIBRARIES
    FFLOW_INCLUDE_DIRS)

mark_as_advanced(FFLOW_INCLUDE_DIR FFLOW_LIBRARY)
