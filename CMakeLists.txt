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

project(icarus_signal_processing VERSION 09.56.00 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

list(APPEND CMAKE_MODULE_PATH
#     ${CMAKE_CURRENT_SOURCE_DIR}/Modules
)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic
  EXTRA_CXX_FLAGS -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( fftw )

# macros for dictionary and simple_plugin
include(CetMake)

# add cet_find_library commands here when needed

##cet_find_library(FFTW_LIBRARY NAMES fftw3f fftw3 fftw3q fftw3l fftw3_threads fftw3f_threads fftw3q_threads PATHS $ENV{FFTW_LIBRARY} )
cet_find_library(ICARUS_FFTW3_LIBRARY  NAMES fftw3  PATHS ENV FFTW_LIBRARY NO_DEFAULT_PATH )
cet_find_library(ICARUS_FFTW3F_LIBRARY NAMES fftw3f PATHS ENV FFTW_LIBRARY NO_DEFAULT_PATH )
cet_find_library(ICARUS_FFTW3L_LIBRARY NAMES fftw3l PATHS ENV FFTW_LIBRARY NO_DEFAULT_PATH )
cet_find_library(ICARUS_FFTW3Q_LIBRARY NAMES fftw3q PATHS ENV FFTW_LIBRARY NO_DEFAULT_PATH )
set(ICARUS_FFTW_LIBRARIES ICARUS_FFTW3_LIBRARY ICARUS_FFTW3F_LIBRARY ICARUS_FFTW3L_LIBRARY ICARUS_FFTW3Q_LIBRARY)
include_directories( $ENV{EIGEN_INC} )

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(icarus_signal_processing)

# tests
add_subdirectory(test)

# ups - table and config files
# must be AFTER all other subdirectories
add_subdirectory(ups)

# packaging utility - must be last
include(UseCPack)
