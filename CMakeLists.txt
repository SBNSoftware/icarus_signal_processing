# ======================================================================
#  icarus_signal_processing main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/icarus_signal_processing/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================
cmake_minimum_required(VERSION 3.19 FATAL_ERROR)


find_package(cetmodules 2.25.05 REQUIRED)

project(icarus_signal_processing VERSION 09.60.01 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic
  EXTRA_CXX_FLAGS -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

# these are minimum required versions, not the actual product versions
find_package( PkgConfig REQUIRED)
find_package( ROOT )
# macros for dictionary and simple_plugin
include(CetMake)

# add cet_find_library commands here when needed

pkg_search_module(FFTW3 REQUIRED IMPORTED_TARGET GLOBAL fftw3f fftw3 fftw3q fftw3l fftw3_threads fftw3f_threeads fftw3q_threads)
set(ICARUS_FFTW_LIBRARIES PkgConfig::FFTW3)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(icarus_signal_processing)

# tests
add_subdirectory(test)

# ups - table and config files
# must be AFTER all other subdirectories

# packaging utility - must be last
cet_cmake_config()
