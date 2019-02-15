cmake_minimum_required(VERSION 2.8.6)
project(casadi C CXX)

set(CASADI_MAJOR_VERSION 3)
set(CASADI_MINOR_VERSION 4)
set(CASADI_PATCH_VERSION 5)
set(CASADI_IS_RELEASE 1)

set(CASADI_VERSION ${CASADI_MAJOR_VERSION}.${CASADI_MINOR_VERSION}.${CASADI_PATCH_VERSION})
if(NOT CASADI_IS_RELEASE)
  set(CASADI_VERSION "${CASADI_VERSION}+")
endif()
set(CASADI_VERSION "${CASADI_VERSION}")
set(PACKAGE_NAME "CasADi")
set(PACKAGE_VERSION ${CASADI_VERSION})
set(PACKAGE_SO_VERSION "3.5")
set(PACKAGE_DESCRIPTION "CasADi is a symbolic framework for automatic differentation and numeric optimization")
set(PACKAGE_AUTHOR "Joel Andersson and Joris Gillis")
set(PACKAGE_MAINTAINER "Joel Andersson and Joris Gillis")
set(PACKAGE_URL "http://casadi.org")

############################################################
####################### Policies ###########################
############################################################
# Enforce cmake2.4+ style behaviour for
# mixed filename/full path linker directives
cmake_policy(SET CMP0003 NEW)

# cmake_policy propagation rules
if(POLICY CMP0011)
cmake_policy(SET CMP0011 OLD)
endif()

# cmake_policy avoid escaping in add_definitions
cmake_policy(SET CMP0005 OLD)

# cmake_policy use MACOSX_RPATH=ON by default
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

# cmake_policy how to interpret string arguments
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

#############################################################
########################## Version ##########################
#############################################################
# For systems with git installed, find out revsion and description.
execute_process(COMMAND git rev-parse HEAD
  OUTPUT_VARIABLE git_revision
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
execute_process(COMMAND git describe --first-parent HEAD
  OUTPUT_VARIABLE git_describe
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
string(REPLACE "+-" "+" git_describe "${git_describe}")
string(REPLACE "-g" "." git_describe "${git_describe}")
set(PACKAGE_VERSION_FULL "${git_describe}")

# overwrite version numbers for official release
if((PACKAGE_VERSION_FULL STREQUAL "") OR (NOT PACKAGE_VERSION MATCHES "\\+"))
  set(PACKAGE_VERSION_FULL "${PACKAGE_VERSION}")
endif()

# We use the FeatureSummary module to get a more readable CMake output
include(FeatureSummary)
if(NOT COMMAND add_feature_info)
  function(add_feature_info _name _enabled _description)
    message(STATUS "Detected: ${_name} ${_description}")
  endfunction()
  function(feature_summary)
  endfunction()
endif()

# Default build type is Release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

########################################################################
############################# Fortran ##################################
########################################################################
option(FORTRAN_REQUIRED "Make fortran a required dependency" OFF)

# http://public.kitware.com/Bug/print_bug_page.php?bug_id=9220
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(language_support)
workaround_9220(Fortran Fortran_language_works)
if(FORTRAN_REQUIRED)
  enable_language(Fortran)
  set(Fortran_language_works ON)
else()
  if(Fortran_language_works)
    enable_language(Fortran OPTIONAL)
  endif()
endif()

########################################################################
############################# C++11 ####################################
########################################################################
# Check if the compiler supports C++11
foreach(CXX11FLAG_TEST "" "-std=c++11" "-std=c++0x")
  if(NOT USE_CXX11)
    # Try to compile test with the proposed flag
    try_compile(USE_CXX11
      ${CMAKE_BINARY_DIR}
      ${CMAKE_MODULE_PATH}/check_cxx11.cpp
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG_TEST}
      OUTPUT_VARIABLE OUTPUT)

    # Save flag, if successful
    if(USE_CXX11)
      set(CXX11FLAG ${CXX11FLAG_TEST} CACHE INTERNAL "C++11 compile flag")
    endif()
  endif()
endforeach()

# If found, set flag
if(USE_CXX11)
  if(CXX11FLAG)
    message(STATUS "Flag needed for enabling C++11 features: ${CXX11FLAG}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11FLAG}")
  else()
    message(STATUS "No flag needed for enabling C++11 features.")
  endif()
  add_definitions(-DUSE_CXX11)
else()
  message(FATAL_ERROR "Cannot enable C++11 features.")
endif()

if(CXX11FLAG)
  try_compile(HAS_COPYSIGN
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_copysign.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_COPYSIGN)
    add_definitions(-DHAS_COPYSIGN)
  endif()


  try_compile(HAS_ERF
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_erf.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_ERF)
    add_definitions(-DHAS_ERF)
  endif()

  try_compile(HAS_ERFINV
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_erfinv.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_ERFINV)
    add_definitions(-DHAS_ERFINV)
  endif()
endif()

# Check if mkstemps is available
if(WIN32)
  add_definitions(-DHAVE_SIMPLE_MKSTEMPS)
else()
  file(WRITE ${CMAKE_BINARY_DIR}/mkstemps_test.cxx
      "#include <cstdlib>\n"
      "#include <unistd.h>\n"
      "int main() {\n"
      "  char tmp[] = \"mkstempstestXXXXXX.so\";\n"
      "  int flag = mkstemps(tmp, 3);\n"
      "  return flag;\n"
      "}\n")
  try_compile(HAVE_MKSTEMPS ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/mkstemps_test.cxx
              OUTPUT_VARIABLE HAVE_MKSTEMPS_OUTPUT)
  file(WRITE ${CMAKE_BINARY_DIR}/simple_mkstemps_test.cxx
      "#include <sys/stat.h>\n"
      "#include <fcntl.h>\n"
      "#include <errno.h>\n"
      "int main() {\n"
      "  int fd = open(\"test.c\", O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);\n"
      "}\n")
  try_compile(HAVE_SIMPLE_MKSTEMPS ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/simple_mkstemps_test.cxx
              OUTPUT_VARIABLE HAVE_SIMPLE_MKSTEMPS_OUTPUT)
  if(HAVE_MKSTEMPS)
    add_definitions(-DHAVE_MKSTEMPS)
  else()
    if(HAVE_SIMPLE_MKSTEMPS)
      add_definitions(-DHAVE_SIMPLE_MKSTEMPS)
    else()
      message(WARNING "mkstemp is not available, falling back on tmpnam. "
                      "Test compilation returned the following:\n"
                      "${HAVE_MKSTEMPS_OUTPUT}"
                      "Test compilation returned the following:\n"
                      "${HAVE_SIMPLE_MKSTEMPS_OUTPUT}")
   endif()
  endif()
endif()

############################################################
################# WITH_COMMON option #######################
############################################################
option(WITH_COMMON "Build common, stable interfaces, if found" OFF)
if(WITH_COMMON)
  # Look for LAPACK
  find_package(LAPACK QUIET)
  set(WITH_LAPACK_DEF ${LAPACK_FOUND})
  # qpOASES and blockSQP built if LAPACK available
  set(WITH_QPOASES_DEF ${LAPACK_FOUND})
  set(WITH_BLOCKSQP_DEF ${LAPACK_FOUND})
  # Look for IPOPT
  find_package(IPOPT QUIET)
  set(WITH_IPOPT_DEF ${IPOPT_FOUND})
else()
  # Disabled by default
  set(WITH_LAPACK_DEF OFF)
  set(WITH_QPOASES_DEF OFF)
  set(WITH_BLOCKSQP_DEF OFF)
  set(WITH_IPOPT_DEF OFF)
endif()


############################################################
################# Misc options #############################
############################################################
option(INSTALL_INTERNAL_HEADERS "Install all headers" OFF)
option(ENABLE_EXPORT_ALL "Export all symbols to a shared library" OFF)
option(ENABLE_STATIC "Compile static libraries (.a on Linux)" OFF)
option(ENABLE_SHARED "Compile shared libraries (.so on Linux)" ON)
option(SWIG_EXPORT "Export SWIG" OFF)
option(SWIG_IMPORT "Export SWIG" OFF)

# enabling OpenMP support if requested
option(WITH_OPENMP "Compile with parallelization support" OFF)
if(WITH_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DWITH_OPENMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DWITH_OPENMP")
  endif()
endif()

# enabling Threads support if requested
option(WITH_THREAD "Compile with parallelization support" OFF)
option(WITH_THREAD_MINGW "For MINGW, use mingw-std-threads" OFF)
if(WITH_THREAD)
  add_definitions(-DCASADI_WITH_THREAD)
endif()
if(MINGW AND WITH_THREAD_MINGW)
  add_definitions(-DCASADI_WITH_THREAD_MINGW)
else()
  if(WITH_THREAD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()
endif()

# Enable: RTLD_DEEPBIND
option(WITH_DEEPBIND "Load plugins with RTLD_DEEPBIND (can be used to resolve conflicting libraries in e.g. MATLAB)" ON)
if(WITH_DEEPBIND)
  add_definitions(-DWITH_DEEPBIND)
endif()

# Suppress warnings for refcounting failures (at program closure)
option(WITH_REFCOUNT_WARNINGS "Issue warnings upon reference counting failure" OFF)
if(WITH_REFCOUNT_WARNINGS)
  add_definitions(-DWITH_REFCOUNT_WARNINGS)
endif()

# Have an so version?
option(WITH_SO_VERSION "Use an so version for the library (version suffix) when applicable" ON)

set(WITH_COPYSIGN_UNDEF_DEFAULT OFF)
if(WITH_PYTHON AND MSVC)
  set(WITH_COPYSIGN_UNDEF_DEFAULT OFF)
endif()
option(WITH_COPYSIGN_UNDEF "On Windows, in Python.h #defines copysign to be _copysign. This option inserts an #undef" WITH_COPYSIGN_UNDEF_DEFAULT)

############################################################
################# Installation location ####################
############################################################
# Flat hierarchy if self-contained and on Windows
option(WITH_SELFCONTAINED "Make the install directory self-contained" OFF)
if(WITH_SELFCONTAINED OR (WIN32 AND NOT CYGWIN))
  set(DEF_LIB_PREFIX casadi)
  set(DEF_BIN_PREFIX casadi)
  set(DEF_CMAKE_PREFIX casadi/cmake)
  set(DEF_PYTHON_PREFIX) # Note: Files installed in a casadi subdirectory
  set(DEF_MATLAB_PREFIX casadi)
  set(DEF_INCLUDE_PREFIX casadi/include)
else()
  set(DEF_LIB_PREFIX lib)
  set(DEF_BIN_PREFIX bin)
  set(DEF_CMAKE_PREFIX lib/cmake/casadi)
  set(DEF_PYTHON_PREFIX python)
  if (WITH_PYTHON)
    set(PYTHON_VERSION_REQ "")
    if (WITH_PYTHON3)
      set(PYTHON_VERSION_REQ "3")
    endif()
    find_package(PythonInterp ${PYTHON_VERSION_REQ} REQUIRED)
    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; from distutils.sysconfig import get_python_lib; print(get_python_lib())"
      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(DEF_PYTHON_PREFIX ${PYTHON_SITE_PACKAGES})
  endif()
  set(DEF_MATLAB_PREFIX matlab)
  set(DEF_INCLUDE_PREFIX include)
endif()

# Installation directories are cache entries that can be modified
set(LIB_PREFIX ${DEF_LIB_PREFIX} CACHE PATH "Installation directory for libraries")
set(BIN_PREFIX ${DEF_BIN_PREFIX} CACHE PATH "Installation directory for executables")
set(INCLUDE_PREFIX ${DEF_INCLUDE_PREFIX} CACHE PATH "Installation directory for header files")
set(CMAKE_PREFIX ${DEF_CMAKE_PREFIX} CACHE PATH "Installation directory for CMake files")
set(PYTHON_PREFIX ${DEF_PYTHON_PREFIX} CACHE PATH "Installation directory for Python files")
set(MATLAB_PREFIX ${DEF_MATLAB_PREFIX} CACHE PATH "Installation directory for MATLAB/Octave files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE PYTHON MATLAB)
  set(var ${p}_PREFIX)
  if(${var})
    if(NOT IS_ABSOLUTE "${${var}}")
      # Turn relative path into absolute
      set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif()
  else()
    # Flat installation hierarchy
    set(${var} "${CMAKE_INSTALL_PREFIX}")
  endif()
endforeach()

# Set rpath
if(WITH_SELFCONTAINED)
  if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/")
  else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
  endif()
else()
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${LIB_PREFIX}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${LIB_PREFIX}")
  endif()
endif()


################################################################################
############################### compiler warnings ##############################
################################################################################
# add compiler warnings
option(WITH_EXTRA_WARNINGS "Add extra warnings (like -Wall -Wextra)" OFF)
option(WITH_WERROR "Add -Werror" OFF)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-attributes")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-maybe-uninitialized -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-mismatched-tags -Wno-tautological-compare")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-maybe-uninitialized -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # I don't know the MSVC warnings, so do nothing for now
endif()
#######################################################################
############################### Extra runtime checks ##################
#######################################################################
option(WITH_EXTRA_CHECKS "Compile in extra runtime checks (for developers)" OFF)
if(WITH_EXTRA_CHECKS)
   add_definitions(-DWITH_EXTRA_CHECKS)
endif()


#######################################################################
############################# -fPIC ###################################
#######################################################################
# the following commands are needed to fix a problem with the libraries
# for linux 64 bits
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
  message(STATUS "x86_64 architecture detected - setting flag -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
endif()

if(MINGW)
  # Circumventing a bug in MinGW g++ v4.7.2, evoked by 752fa89355ffa
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-ipa-cp-clone")
endif()

# Tentative fix for https://github.com/casadi/casadi/issues/1665
if(MINGW OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel"))
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fkeep-inline-functions")
endif()

# Plugin version
add_definitions(-DCASADI_VERSION=31)

# Necessary for Visual C++
add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)

# Dynamic loading (needed for loading plugins and external functions)
option(WITH_DL "Enable dynamic loading of functions" ON)
if(WITH_DL)
  add_definitions(-DWITH_DL)
endif()
add_feature_info(dynamic-loading WITH_DL "Compile with support for dynamic loading of generated functions (needed for ExternalFunction)")

# Include support for deprecated features (to be removed in the next release)
option(WITH_DEPRECATED_FEATURES "Compile with syntax that is scheduled to be deprecated" ON)
if (WITH_DEPRECATED_FEATURES)
  add_definitions(-DWITH_DEPRECATED_FEATURES)
endif()

include_directories(.)
include_directories(${PROJECT_BINARY_DIR})

######################################################
##################### paths ##########################
######################################################
if(WIN32)
  # On Windows, build everything in the same directory
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
else()
  # On non-Windows, use separate directories for build
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
endif()

#####################################################
################# include the rest ##################
#####################################################
add_subdirectory(external_packages)
add_subdirectory(casadi)
add_subdirectory(experimental EXCLUDE_FROM_ALL)
add_subdirectory(misc)

option(WITH_EXAMPLES "Build examples" ON)
if(WITH_EXAMPLES)
  add_subdirectory(docs/examples)
  add_subdirectory(docs/api/examples/ctemplate)
endif()


########################################################
########### Generate config files ######################
########################################################
# Add all targets to the build-tree export set
export(TARGETS casadi
  FILE "${PROJECT_BINARY_DIR}/casadi-targets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE CasADi)

# Create the casadi-config.cmake and casadi-config-version.cmake files
file(RELATIVE_PATH REL_INCLUDE_PREFIX "${CMAKE_PREFIX}" "${INCLUDE_PREFIX}")
file(RELATIVE_PATH REL_LIB_PREFIX "${CMAKE_PREFIX}" "${LIB_PREFIX}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
set(CONF_LIBRARY_DIRS "${LIBRARY_OUTPUT_PATH}")
configure_file(casadi-config.cmake.in "${PROJECT_BINARY_DIR}/casadi-config.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${CASADI_CMAKE_DIR}/${REL_INCLUDE_PREFIX}")
set(CONF_LIBRARY_DIRS "\${CASADI_CMAKE_DIR}/${REL_LIB_PREFIX}")
configure_file(casadi-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/casadi-config.cmake" @ONLY)
# ... for both
configure_file(casadi-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/casadi-config-version.cmake" @ONLY)

# Install the casadi-config.cmake and casadi-config-version.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/casadi-config.cmake"
  "${PROJECT_BINARY_DIR}/casadi-config-version.cmake"
  DESTINATION "${CMAKE_PREFIX}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT casadi-targets DESTINATION "${CMAKE_PREFIX}" COMPONENT dev)

set(CPACK_PACKAGE_CONTACT "casadi-users@googlegroups.com")
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION_FULL})
include(CPack)
