# Install script for directory: /workspace

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/fflow" TYPE FILE FILES
    "/workspace/include/fflow/config.hh"
    "/workspace/include/fflow/common.hh"
    "/workspace/include/fflow/small_vector.hh"
    "/workspace/include/fflow/refcounted_ptr.hh"
    "/workspace/include/fflow/shared_array.hh"
    "/workspace/include/fflow/debug.hh"
    "/workspace/include/fflow/function_cache.hh"
    "/workspace/include/fflow/gcd.hh"
    "/workspace/include/fflow/integer_math.hh"
    "/workspace/include/fflow/matrix.hh"
    "/workspace/include/fflow/multivariate_reconstruction.hh"
    "/workspace/include/fflow/multivariate_reconstruction_details.hh"
    "/workspace/include/fflow/polynomial.hh"
    "/workspace/include/fflow/primes.hh"
    "/workspace/include/fflow/rational_function.hh"
    "/workspace/include/fflow/univariate_reconstruction.hh"
    "/workspace/include/fflow/ratfun_parser.hh"
    "/workspace/include/fflow/mp_common.hh"
    "/workspace/include/fflow/mp_functions.hh"
    "/workspace/include/fflow/mp_gcd.hh"
    "/workspace/include/fflow/mp_multivariate_reconstruction.hh"
    "/workspace/include/fflow/algorithm.hh"
    "/workspace/include/fflow/alg_linear_solver.hh"
    "/workspace/include/fflow/alg_linear_fit.hh"
    "/workspace/include/fflow/alg_reconstruction.hh"
    "/workspace/include/fflow/alg_mp_reconstruction.hh"
    "/workspace/include/fflow/format.h"
    "/workspace/include/fflow/ostream.h"
    "/workspace/include/fflow/json.hh"
    "/workspace/include/fflow/thread_pool.hh"
    "/workspace/include/fflow/graph.hh"
    "/workspace/include/fflow/subgraph.hh"
    "/workspace/include/fflow/alg_lists.hh"
    "/workspace/include/fflow/alg_functions.hh"
    "/workspace/include/fflow/analytic_solver.hh"
    "/workspace/include/fflow/numeric_solver.hh"
    "/workspace/include/fflow/node_solver.hh"
    "/workspace/include/fflow/analytic_fit.hh"
    "/workspace/include/fflow/subgraph_fit.hh"
    "/workspace/include/fflow/numeric_fit.hh"
    "/workspace/include/fflow/alg_laurent.hh"
    "/workspace/include/fflow/subgraph_reconstruct.hh"
    "/workspace/include/fflow/cached_subgraph.hh"
    "/workspace/include/fflow/capi.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libfflow.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libfflow.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libfflow.so"
         RPATH "/usr/local/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/workspace/fresh_build/libfflow.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libfflow.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libfflow.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libfflow.so"
         OLD_RPATH "::::::::::::::"
         NEW_RPATH "/usr/local/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libfflow.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/workspace/fresh_build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
