# This file defines the macro atlas4py_add_atlas, which is responsible for ensuring that the atlas library
# is available for building the Python bindings.
# It does so by first trying to find an existing installation of atlas using CMake's find_package.
# If it cannot find it, it uses FetchContent to download and build atlas and its dependencies (eckit and ecbuild) from source.
# The macro is called in the main CMakeLists.txt file of the project.

function( atlas4py_parse_version version_str )
  set( options )
  set( single_value_args PREFIX )
  set( multi_value_args )

  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}" ${ARGN} )

  if( NOT _PAR_PREFIX )
    set( prefix "_" )
  else()
    set( prefix ${_PAR_PREFIX} )
  endif()

  ## Parse version_str
  set( ${prefix}_VERSION_STR "${version_str}" )
  string( REGEX REPLACE "^((([0-9]+)\\.)+([0-9]+)).*" "\\1" ${prefix}_VERSION "${version_str}" )
  string( LENGTH "${${prefix}_VERSION}" ver_len )
  string( SUBSTRING "${version_str}" ${ver_len} -1 ${prefix}_VERSION_SUFFIX )
  string( REPLACE "." " " _version_list ${${prefix}_VERSION} ) # dots to spaces
  separate_arguments( _version_list )
  list (LENGTH _version_list _len)
  if( ${_len} GREATER 0 )
    list( GET _version_list 0 ${prefix}_VERSION_MAJOR )
  endif()
  if( ${_len} GREATER 1 )
    list( GET _version_list 1 ${prefix}_VERSION_MINOR )
  endif()
  if( ${_len} GREATER 2 )
    list( GET _version_list 2 ${prefix}_VERSION_PATCH )
  endif()
  if( ${_len} GREATER 3 )
    list( GET _version_list 3 ${prefix}_VERSION_TWEAK )
  endif()

  math(EXPR ${prefix}_VERSION_INT "${${prefix}_VERSION_MAJOR} * 10000 + ${${prefix}_VERSION_MINOR} * 100 + ${${prefix}_VERSION_PATCH}")

  ## Export variables to parent scope
  list( APPEND export_variables_parent_scope
    ${prefix}_VERSION_STR
    ${prefix}_VERSION
    ${prefix}_VERSION_MAJOR
    ${prefix}_VERSION_MINOR
    ${prefix}_VERSION_PATCH
    ${prefix}_VERSION_TWEAK
    ${prefix}_VERSION_SUFFIX
    ${prefix}_VERSION_INT
  )
  foreach( _var ${export_variables_parent_scope} )
    if( DEFINED ${_var} )
      set( ${_var} ${${_var}} PARENT_SCOPE )
    endif()
  endforeach()

endfunction()

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
        find_package(atlas CONFIG REQUIRED)
    endif()

    include(cmake/atlas4py_parse_version.cmake)
    atlas4py_parse_version(${atlas_VERSION_STR} PREFIX ATLAS4PY_ATLAS)
endmacro()
