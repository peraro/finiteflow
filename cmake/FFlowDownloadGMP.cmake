# This is used when -DFFLOW_DOWNLOAD_GMP=1

if(CMAKE_VERSION VERSION_LESS 3.14)
    message(FATAL_ERROR "The option FFLOW_DOWNLOAD_GMP required CMake>=3.14")
endif()

include(ExternalProject)
include(ProcessorCount)

if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()
set(FFLOW_GMP_VERSION 6.3.0)
set(FFLOW_GMP_URLS
    "https://ftp.gnu.org/gnu/gmp/gmp-${FFLOW_GMP_VERSION}.tar.xz"
    "https://ftpmirror.gnu.org/gmp/gmp-${FFLOW_GMP_VERSION}.tar.xz"
    "https://gmplib.org/download/gmp/gmp-${FFLOW_GMP_VERSION}.tar.xz"
)
set(FFLOW_GMP_SHA256
    "a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898")

ProcessorCount(_fflow_processor_count)
if(_fflow_processor_count EQUAL 0)
    set(_fflow_processor_count 1)
endif()

set(_fflow_gmp_lib "${PROJECT_BINARY_DIR}/gmp_install/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX}")

ExternalProject_Add(
    fflow_external_gmp
    URL ${FFLOW_GMP_URLS}
    PREFIX "${PROJECT_BINARY_DIR}/gmp"
    INSTALL_DIR "${PROJECT_BINARY_DIR}/gmp_install"
    CONFIGURE_COMMAND
    <SOURCE_DIR>/configure
    --prefix=<INSTALL_DIR>
    --libdir=<INSTALL_DIR>/lib
    --disable-shared --enable-static --with-pic --disable-cxx
    "CC=cc -std=gnu99"
    "MAKEINFO=true" # prevents running the real makeinfo
    BUILD_COMMAND make -j${_fflow_processor_count}
    INSTALL_COMMAND make install
    BUILD_BYPRODUCTS ${_fflow_gmp_lib}
)

file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/gmp_install/include")
add_library(fflow_downloaded_gmp STATIC IMPORTED)
set_target_properties(fflow_downloaded_gmp
    PROPERTIES IMPORTED_LOCATION "${_fflow_gmp_lib}")
add_dependencies(fflow_downloaded_gmp fflow_external_gmp)

# These are needed by CMakeLists.txt
set(GMP_INCLUDE_DIRS "${PROJECT_BINARY_DIR}/gmp_install/include")
set(GMP_TARGETS fflow_downloaded_gmp)
