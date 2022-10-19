add_definitions(-DEIGEN_FFTW_DEFAULT)
include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})
cet_make(  ### MIGRATE-ACTION-RECOMMENDED: use cet_make_library(), build_dictionary(), basic_plugin() with explicit source lists
          LIBRARIES 
                          ${ICARUS_FFTW_LIBRARIES}
                          ROOT::Geom
                          ROOT::XMLIO
                          ROOT::Gdml
                          ROOT::FFTW
                          ROOT::Core
       )

install_headers()
install_fhicl()
install_source()
