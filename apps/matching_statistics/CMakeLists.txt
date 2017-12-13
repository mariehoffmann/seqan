# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /matching_statistics/CMakeLists.txt
#
# CMakeLists.txt file for matching_statistics.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_matching_statistics CXX)
message (STATUS "Configuring apps/matching_statistics")
# TODO: which version?
set (SEQAN_APP_VERSION "2.3.0")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP)
    find_package (ZLIB)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

if (NOT ZLIB_FOUND)
  message (STATUS "  ZLIB not found: matching_statistics will not be built.")
  return ()
endif (NOT ZLIB_FOUND)

# Enable C++11 in GCC and/or CLANG.
#if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
#    add_definitions ("-std=c++0x")
#endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# We define a library for the reusable parts of Mason.
add_library (ms_lib
             unidirectional_matching_statistics.h)

target_link_libraries (ms_lib ${SEQAN_LIBRARIES})

# Tests for the library functionalities.
add_executable (ms_tests matching_statistics_tests.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (ms_tests ${SEQAN_LIBRARIES} ms_lib)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include mason executables in CTD structure.
set (SEQAN_CTD_EXECUTABLES
    ${SEQAN_CTD_EXECUTABLES}
    ms_tests
    CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (matching_statistics)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install Mason binaries in ${PREFIX}/bin directory
install (TARGETS ms_tests DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/matching_statistics for SeqAn release builds.
install (FILES LICENSE
               README
               README.ms_tests
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/adeno_virus.fa
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (matching_statistics)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:matching_statistics")
  set (CPACK_PACKAGE_NAME "matching_statistics")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Matching statistics package for sequence comparison.")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Marie Hoffmann <marie.hoffmann@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Marie Hoffmann <marie.hoffmann@fu-berlin.de>")

  seqan_configure_cpack_app (matching_statistics "matching_statistics")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:matching_statistics")
