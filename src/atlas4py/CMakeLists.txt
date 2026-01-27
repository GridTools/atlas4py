cmake_minimum_required(VERSION 3.25.0)

if (NOT SKBUILD)
    set (PROJECT_ROOT_DIR "../..")
    cmake_path(ABSOLUTE_PATH PROJECT_ROOT_DIR BASE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR} NORMALIZE OUTPUT_VARIABLE PROJECT_ROOT_DIR)

    message(FATAL_ERROR "\

    This CMake file is meant to be executed using 'scikit-build'. Running
    it directly will almost certainly not produce the desired result. If
    you are a user trying to install this package, please use the command
    below, which will install all necessary build dependencies, compile
    the package in an isolated environment, and then install it.
    =====================================================================
    $ pip install ${PROJECT_ROOT_DIR}
    =====================================================================
    For development purposes, it is usually much more efficient to
    install the build dependencies in your environment once and use the
    following command that avoids a costly creation of a new virtual
    environment at every compilation:
    =====================================================================
    $ pip install pybind11 scikit-build-core
    $ pip install --no-build-isolation -ve ${PROJECT_ROOT_DIR}
    =====================================================================
    You may optionally add -Ceditable.rebuild=true to auto-rebuild when
    the package is imported. Otherwise, you need to re-run the above
    after editing C++ files.")
endif()

set( PROJECT_NAME ${SKBUILD_PROJECT_NAME} )
set( PROJECT_VERSION ${SKBUILD_PROJECT_VERSION} )
set( PROJECT_VERSION_FULL ${SKBUILD_PROJECT_VERSION_FULL} )
string( REPLACE "${PROJECT_VERSION}.dev" "" PROJECT_VERSION_DEV "${PROJECT_VERSION_FULL}" )

message( STATUS "[${PROJECT_NAME}] (${PROJECT_VERSION_FULL})" )
project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})

# Policies
cmake_policy(SET CMP0074 NEW)
if (POLICY CMP0148)
    cmake_policy(SET CMP0148 NEW)
endif()
if (POLICY CMP0168)
    cmake_policy(SET CMP0168 NEW)
endif()

set(CMAKE_CXX_STANDARD 17)

include(cmake/atlas4py_add_atlas.cmake)
atlas4py_add_atlas()

### Find pybind11

message( STATUS "${PROJECT_NAME}: find_package(pybind11 CONFIG)..." )
find_package(pybind11 CONFIG)
if (NOT pybind11_FOUND)
    message( FATAL_ERROR "pybind11 not found. Please install pybind11 or use pip to install this package." )
endif()

### RPATH handling

include(GNUInstallDirs) # defines CMAKE_INSTALL_LIBDIR
if(APPLE)
    set(rpath_origin_install_libdir "@loader_path/${CMAKE_INSTALL_LIBDIR}")
else()
    set(rpath_origin_install_libdir "$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
endif()
set( CMAKE_INSTALL_RPATH "${rpath_origin_install_libdir}" ) # A semicolon-separated list specifying the RPATH to use in installed targets
set( CMAKE_SKIP_BUILD_RPATH              FALSE ) # use RPATHs for the build tree
set( CMAKE_BUILD_WITH_INSTALL_RPATH      TRUE  ) # build with *relative* rpaths by default
if (NOT build_atlas)
  set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE  ) # add the automatic parts to RPATH which point to dirs outside build tree, i.e. external atlas
endif()

### Python bindings module atlas4py

pybind11_add_module(_atlas4py _atlas4py.cpp)
target_link_libraries(_atlas4py PUBLIC atlas)
target_compile_definitions(_atlas4py PRIVATE ATLAS4PY_VERSION_STRING=${PROJECT_VERSION_FULL})

### Installation

install(TARGETS _atlas4py DESTINATION .)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION .
    FILES_MATCHING PATTERN "*.py"
)
