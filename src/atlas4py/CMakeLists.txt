cmake_minimum_required(VERSION 3.25.0)

project(atlas4py LANGUAGES CXX)

cmake_policy(SET CMP0074 NEW)
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

find_package(atlas QUIET PATHS $ENV{ATLAS_INSTALL_DIR})
if( atlas_FOUND )
  message( "Found atlas: ${atlas_DIR} (found version \"${atlas_VERSION}\")" )
endif()

if(NOT atlas_FOUND)
    find_package(ecbuild)
    if(NOT ecbuild_FOUND)
        FetchContent_Declare(
            ecbuild
            GIT_REPOSITORY https://github.com/ecmwf/ecbuild.git
            GIT_TAG        ${ATLAS4PY_ECBUILD_VERSION}
        )
        FetchContent_MakeAvailable(ecbuild)
    endif()
    find_package(eckit)
    if( eckit_FOUND )
      message( "Found eckit: ${eckit_DIR} (found version \"${eckit_VERSION}\")" )
    endif()

    if(NOT eckit_FOUND)
        FetchContent_Declare(
            eckit
            GIT_REPOSITORY https://github.com/ecmwf/eckit.git
            GIT_TAG        ${ATLAS4PY_ECKIT_VERSION}
        )
        FetchContent_MakeAvailable(eckit)
        set(_atlas4py_built_eckit ON)
    endif()

    FetchContent_Declare(
        atlas
        GIT_REPOSITORY https://github.com/ecmwf/atlas.git
        GIT_TAG        ${ATLAS4PY_ATLAS_VERSION}
    )
    set( ENABLE_GRIDTOOLS_STORAGE OFF CACHE BOOL "" FORCE )
    FetchContent_MakeAvailable(atlas)
endif()

find_package(pybind11)

if(APPLE)
    set(rpath_origin_install_libdir "@loader_path/${CMAKE_INSTALL_LIBDIR}")
else()
    set(rpath_origin_install_libdir "$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
endif()

pybind11_add_module(_atlas4py _atlas4py.cpp)
target_link_libraries(_atlas4py PUBLIC atlas)
install(TARGETS _atlas4py DESTINATION .)
set_target_properties(_atlas4py PROPERTIES INSTALL_RPATH "${rpath_origin_install_libdir}")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION .
    FILES_MATCHING PATTERN "*.py"
)
