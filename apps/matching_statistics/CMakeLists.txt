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
set (SEQAN_APP_VERSION "2.0.8")
set(CMAKE_CXX_STANDARD 14)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (Boost COMPONENTS filesystem REQUIRED)
    find_package (OpenMP)
    find_package (ZLIB)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

if (NOT Boost_FOUND OR NOT ZLIB_FOUND OR NOT OPENMP_FOUND)
  message (STATUS "  Boost, ZLIB or OpenMP not found: matching_statistics will not be built.")
  return ()
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (SYSTEM)
include_directories (${Boost_INCLUDE_DIRS})


# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
#" -lboost-filesystem -lboost_system")

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Library for reusable parts of the matching statistics.
add_library (ms_lib
            main.cpp
            unidirectional_matching_statistics.h)

# Add tests.
add_executable (ms_tests matching_statistics_test.cpp)
target_compile_options(ms_tests PUBLIC -std=c++14 -fpermissive -w -Wall)
# maybe only needed when ms_lib needs the seqan libs
set_target_properties(ms_tests PROPERTIES LINKER_LANGUAGE CXX)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (ms_tests LINK_PUBLIC ${SEQAN_LIBRARIES} ${Boost_LIBRARIES} ms_lib)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include mason executables in CTD structure.
set (SEQAN_CTD_EXECUTABLES
    ${SEQAN_CTD_EXECUTABLES}
    ms_tests
    CACHE INTERNAL "")

# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

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
#install (FILES example/adeno_virus.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

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
