# This file defines the macro atlas4py_add_atlas, which is responsible for ensuring that the atlas library
# is available for building the Python bindings.
# It does so by first trying to find an existing installation of atlas using CMake's find_package.
# If it cannot find it, it uses FetchContent to download and build atlas and its dependencies (eckit and ecbuild) from source.
# The macro is called in the main CMakeLists.txt file of the project.

macro(atlas4py_add_atlas)
    if (NOT build_atlas)
        if (DEFINED ENV{ATLAS_INSTALL_DIR})
            set( atlas_ROOT $ENV{ATLAS_INSTALL_DIR} )
        endif()
        find_package(atlas CONFIG QUIET)
        if( atlas_FOUND )
            message( STATUS "Found atlas: ${atlas_DIR} (found version \"${atlas_VERSION}\")" )
            message( STATUS "Found eckit: ${eckit_DIR} (found version \"${eckit_VERSION}\")" )
            message( STATUS "If instead you want ensure to build atlas from source, configure with -Dbuild_atlas=ON")
            if (NOT atlas_VERSION VERSION_EQUAL ATLAS4PY_ATLAS_VERSION)
                message( WARNING "Found atlas version \"${atlas_VERSION}\", but configured version is \"${ATLAS4PY_ATLAS_VERSION}\"" )
            endif()
            if (NOT eckit_VERSION VERSION_EQUAL ATLAS4PY_ECKIT_VERSION)
                message( WARNING "Found eckit version \"${eckit_VERSION}\", but configured version is \"${ATLAS4PY_ECKIT_VERSION}\"" )
            endif()
            set( ignore_variable ${ATLAS4PY_ECBUILD_VERSION} ) # To avoid "unused variable" warning
        else()
            set(build_atlas TRUE CACHE INTERNAL "build atlas")
        endif()
    endif()
    if(build_atlas)
        include(FetchContent)

        ### Download dependencies
        set ( ecbuild_ROOT ${CMAKE_BINARY_DIR}/_deps/ecbuild )
        message( STATUS "Downloading ecbuild version \"${ATLAS4PY_ECBUILD_VERSION}\" to ${ecbuild_ROOT}" )
        message( STATUS "Downloading and building eckit version \"${ATLAS4PY_ECKIT_VERSION}\"" )
        message( STATUS "Downloading and building atlas version \"${ATLAS4PY_ATLAS_VERSION}\"" )
        FetchContent_Populate(
            ecbuild
            GIT_REPOSITORY https://github.com/ecmwf/ecbuild.git
            GIT_TAG        ${ATLAS4PY_ECBUILD_VERSION}
            SOURCE_DIR     ${ecbuild_ROOT}
            QUIET
        )
        FetchContent_Declare(
            eckit
            GIT_REPOSITORY https://github.com/ecmwf/eckit.git
            GIT_TAG        ${ATLAS4PY_ECKIT_VERSION}
        )
        FetchContent_Declare(
            atlas
            GIT_REPOSITORY https://github.com/ecmwf/atlas.git
            GIT_TAG        ${ATLAS4PY_ATLAS_VERSION}
        )

        # Disable unused features for faster compilation
        set(ECKIT_ENABLE_TESTS     OFF)
        set(ECKIT_ENABLE_DOCS      OFF)
        set(ECKIT_ENABLE_PKGCONFIG OFF)
        set(ECKIT_ENABLE_ECKIT_GEO OFF)
        set(ECKIT_ENABLE_ECKIT_SQL OFF)
        set(ECKIT_ENABLE_ECKIT_CMD OFF)

        set(ATLAS_ENABLE_TESTS     OFF)
        set(ATLAS_ENABLE_DOCS      OFF)
        set(ATLAS_ENABLE_PKGCONFIG OFF)
        set(ECKIT_ENABLE_ECKIT_GEO OFF)
        set(ECKIT_ENABLE_ECKIT_SQL OFF)
        set(ATLAS_ENABLE_ECKIT_CMD OFF)

        FetchContent_MakeAvailable(eckit atlas)
    endif()
endmacro()
