cmake_minimum_required(VERSION 3.4)

project(finiteflow CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release
        CACHE STRING "Build type: Debug|Release" FORCE)
endif()

option(BUILD_SHARED_LIBS "build shared libraries" ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 99)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

set(AMFF_ENABLE_EXCEPTIONS OFF CACHE BOOL
    "Do not add compiler flags disabling exceptions")

if(NOT(AMFF_ENABLE_EXCEPTIONS))
    if(MSVC)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/EHs-c->)
    else()
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions>)
    endif()
endif()

set(FFLOW_THREAD_POOL OFF CACHE BOOL
    "Use thread pools for parallel tasks")


# rpaths (see https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH)

if (POLICY CMP0068)
    cmake_policy(SET CMP0068 NEW)
endif()
if(APPLE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")


# external dependencies

find_package(FLINT REQUIRED)
include_directories(SYSTEM ${FLINT_INCLUDE_DIRS})
include_directories(SYSTEM ${FLINT_INCLUDE_DIRS}/flint)
set(LIBS ${LIBS} ${FLINT_TARGETS})

find_package(GMP REQUIRED)
include_directories(SYSTEM ${GMP_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GMP_TARGETS})

find_package(Threads)
set(LIBS ${LIBS} ${CMAKE_THREAD_LIBS_INIT})


# fflow library

include(TestBigEndian)
TEST_BIG_ENDIAN(FFLOW_BIG_ENDIAN)
configure_file(include/fflow/config.hh.in include/fflow/config.hh)

include_directories(
    BEFORE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

set(FFLOW_SOURCES
    src/debug.cc
    src/function_cache.cc
    src/gcd.cc
    src/integer_math.cc
    src/matrix.cc
    src/multivariate_reconstruction.cc
    src/polynomial.cc
    src/primes.cc
    src/rational_function.cc
    src/univariate_reconstruction.cc
    src/algorithm.cc
    src/alg_linear_solver.cc
    src/alg_linear_fit.cc
    src/alg_reconstruction.cc
    src/alg_mp_reconstruction.cc
    src/ratfun_parser.cc

    src/mp_functions.cc
    src/mp_gcd.cc
	src/mp_multivariate_reconstruction.cc

    src/format.cc
    src/ostream.cc
    src/spooky_v2.hh
    src/spooky_v2.cc

    #src/interface.cc
    src/json.cc
    src/graph.cc
    src/subgraph.cc
    src/alg_lists.cc
    src/alg_functions.cc
    src/analytic_solver.cc
    src/analytic_fit.cc
    src/numeric_solver.cc
    src/node_solver.cc
    src/subgraph_fit.cc
    src/numeric_fit.cc
    src/alg_laurent.cc
    src/subgraph_reconstruct.cc
    src/cached_subgraph.cc
    src/capi.cc
    )

set(FFLOW_HEADERS
    include/fflow/config.hh

    include/fflow/common.hh
    include/fflow/small_vector.hh
    include/fflow/refcounted_ptr.hh
    include/fflow/shared_array.hh
    include/fflow/debug.hh
    include/fflow/function_cache.hh
    include/fflow/gcd.hh
    include/fflow/integer_math.hh
    include/fflow/matrix.hh
    include/fflow/multivariate_reconstruction.hh
    include/fflow/multivariate_reconstruction_details.hh
    include/fflow/polynomial.hh
    include/fflow/primes.hh
    include/fflow/rational_function.hh
    include/fflow/univariate_reconstruction.hh
    include/fflow/ratfun_parser.hh

    include/fflow/mp_common.hh
    include/fflow/mp_functions.hh
    include/fflow/mp_gcd.hh
	include/fflow/mp_multivariate_reconstruction.hh

    include/fflow/algorithm.hh
    include/fflow/alg_linear_solver.hh
    include/fflow/alg_linear_fit.hh
    include/fflow/alg_reconstruction.hh
    include/fflow/alg_mp_reconstruction.hh

    include/fflow/format.h
    include/fflow/ostream.h

    #include/fflow/interface.hh
    include/fflow/json.hh
    include/fflow/thread_pool.hh
    include/fflow/graph.hh
    include/fflow/subgraph.hh
    include/fflow/alg_lists.hh
    include/fflow/alg_functions.hh
    include/fflow/analytic_solver.hh
    include/fflow/numeric_solver.hh
    include/fflow/node_solver.hh
    include/fflow/analytic_fit.hh
    include/fflow/subgraph_fit.hh
    include/fflow/numeric_fit.hh
    include/fflow/alg_laurent.hh
    include/fflow/subgraph_reconstruct.hh
    include/fflow/cached_subgraph.hh
    include/fflow/capi.h
    )

add_library(fflow ${FFLOW_SOURCES})
target_link_libraries(fflow ${LIBS})

install(FILES ${FFLOW_HEADERS} DESTINATION include/fflow)
install(TARGETS fflow ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)


# tests

set(TEST_TARGETS
    testjson
    )

set(TEST_CTARGETS
    testcapi
    )

foreach(target ${TEST_TARGETS})
    add_executable(${target} tests/${target}.cc)
    set_target_properties(${target} PROPERTIES EXCLUDE_FROM_ALL TRUE)
    target_link_libraries(${target} fflow ${LIBS})
endforeach(target)

foreach(target ${TEST_CTARGETS})
    add_executable(${target} tests/${target}.c)
    set_target_properties(${target} PROPERTIES EXCLUDE_FROM_ALL TRUE)
    target_link_libraries(${target} fflow ${LIBS})
endforeach(target)

if(NOT TARGET tests)
  add_custom_target(tests)
  add_dependencies(tests ${TEST_TARGETS} ${TEST_CTARGETS})
endif()


# mathematica

set (CMAKE_MODULE_PATH
    "${CMAKE_SOURCE_DIR}/cmake/Mathematica" ${CMAKE_MODULE_PATH})

set(MATHLIBINSTALL "" CACHE PATH
    "Installation path for Mathematica interface")

set(FFLOW_MATHEMATICA_SOURCES
    mathlink/fflowmlink.cc
    )

find_package(Mathematica COMPONENTS MathLink WolframLibrary)

if (Mathematica_FOUND)

    include_directories(${Mathematica_INCLUDE_DIRS})
    Mathematica_ADD_LIBRARY(fflowmlink ${FFLOW_MATHEMATICA_SOURCES})
    target_link_libraries(fflowmlink fflow ${LIBS}
        ${Mathematica_MathLink_LIBRARIES})
    if (Mathematica_MathLink_LINKER_FLAGS)
		set_target_properties(fflowmlink PROPERTIES
            LINK_FLAGS "${Mathematica_MathLink_LINKER_FLAGS}")
	endif()
    Mathematica_ABSOLUTIZE_LIBRARY_DEPENDENCIES(fflowmlink)
    if (MATHLIBINSTALL)
        install(TARGETS fflowmlink
            LIBRARY DESTINATION "${MATHLIBINSTALL}")
    elseif (DEFINED Mathematica_USERBASE_DIR)
        foreach (_systemID ${Mathematica_SYSTEM_IDS})
            install(TARGETS fflowmlink
                LIBRARY DESTINATION "${Mathematica_USERBASE_DIR}/SystemFiles/LibraryResources/${_systemID}")
        endforeach()
    else()
        message(WARNING "\
The Mathematica interface cannot be installed \
because no suitable installation path was found.  \
Please specify a path with the MATHLIBINSTALL variable\
during configuration."
            )
    endif()

else()

    message(WARNING "\
The Mathematica interface cannot be installed \
because no suitable Mathematica installation was found."
        )

endif()


# uninstall
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
