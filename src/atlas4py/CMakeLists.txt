cmake_minimum_required(VERSION ${ATLAS4PY_CMAKE_MINIMUM_REQUIRED_VERSION})

project(atlas4py LANGUAGES CXX)

include(FetchContent)
macro(find_ecbuild)
    if(NOT ecbuild_FOUND)
        FetchContent_Declare(
            ecbuild
            GIT_REPOSITORY https://github.com/ecmwf/ecbuild.git
            GIT_TAG        ${ATLAS4PY_ECBUILD_VERSION}
        )
        FetchContent_MakeAvailable(ecbuild)
    endif()
endmacro()

function(copy_target ct_tgt)
    add_custom_target(copy_${ct_tgt} ALL
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${ct_tgt}> ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        DEPENDS ${ct_tgt}
    )
endfunction()

find_package(atlas QUIET PATHS $ENV{ATLAS_INSTALL_DIR})
if(NOT atlas_FOUND)
    find_package(eckit QUIET)
    if(NOT eckit_FOUND)
        find_ecbuild()
        FetchContent_Declare(
            eckit
            GIT_REPOSITORY https://github.com/ecmwf/eckit.git
            GIT_TAG        ${ATLAS4PY_ECKIT_VERSION}
        )
        FetchContent_MakeAvailable(eckit)

        set(_atlas4py_built_eckit ON)
    endif()

    find_ecbuild()
    FetchContent_Declare(
        atlas
        GIT_REPOSITORY https://github.com/ecmwf/atlas.git
        GIT_TAG        ${ATLAS4PY_ATLAS_VERSION}
    )
    FetchContent_MakeAvailable(atlas)

    copy_target(atlas)

    if(_atlas4py_built_eckit)
        # copy atlas' eckit dependencies
        get_target_property(deps atlas LINK_LIBRARIES)
        get_target_property(deps_interface atlas INTERFACE_LINK_LIBRARIES)
        list(APPEND deps ${deps_interface})
        list(REMOVE_DUPLICATES deps)
        foreach(dep IN LISTS deps)
            if(dep MATCHES "eckit.*")
                copy_target(${dep})
            endif()
        endforeach()
    endif()
endif()

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v${ATLAS4PY_PYBIND11_VERSION}
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(_atlas4py _atlas4py.cpp)
target_compile_features(_atlas4py PUBLIC cxx_std_17)
target_link_libraries(_atlas4py PUBLIC atlas eckit)
