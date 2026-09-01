# This file defines a function to parse version strings and export version components as CMake variables.

##############################################################################
#.rst:
#
# atlas4py_parse_version
# =====================
#
# Parse version string of the form "<major>[.<minor>[.<patch>[.<tweak>]]][<suffix>]" ::
#
#   atlas4py_parse_version( <version_str> [ PREFIX <prefix> ] )
#
# Options
# -------
#
# PREFIX : optional
#   string to be prefixed to all defined variables. If not given,
#   the value "_" will be used.
#
# Notes
# -----
#
# Following variables if possible:
#
#      <prefix>_VERSION_STR     = <major>[.<minor>[.<patch>[.<tweak>]]][<suffix>]
#      <prefix>_VERSION         = <major>[.<minor>[.<patch>[.<tweak>]]]
#      <prefix>_VERSION_MAJOR   = <major>
#      <prefix>_VERSION_MINOR   = <minor>
#      <prefix>_VERSION_PATCH   = <patch>
#      <prefix>_VERSION_TWEAK   = <tweak>
#      <prefix>_VERSION_SUFFIX  = <suffix>
#      <prefix>_VERSION_INT     = <major> * 10000 + <minor> * 100 + <patch>
#
# Note, this is taken from ecbuild_parse_version.cmake, with addition of <prefix>_VERSION_INT

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

  if( NOT DEFINED ${prefix}_VERSION_MINOR )
    set( ${prefix}_VERSION_MINOR 0 )
  endif()
  if( NOT DEFINED ${prefix}_VERSION_PATCH )
    set( ${prefix}_VERSION_PATCH 0 )
  endif()
  if( NOT DEFINED ${prefix}_VERSION_TWEAK )
    set( ${prefix}_VERSION_TWEAK 0 )
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
