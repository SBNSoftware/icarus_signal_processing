
FILE(GLOB src_files *.cxx)
cet_make_library(
          SOURCE ${src_files}
          LIBRARIES
                          icarus_signal_processing_Filters
                          ROOT::Geom
                          ROOT::XMLIO
                          ROOT::Gdml
                          ROOT::FFTW
                          ROOT::Core
       )

install_headers()
install_fhicl()
install_source()
