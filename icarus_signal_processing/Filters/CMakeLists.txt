FILE(GLOB src_files *.cxx)

cet_make_library(
          SOURCE ${src_files}
          LIBRARIES
                          FFTW3f::FFTW3f
                          ROOT::Geom
                          ROOT::XMLIO
                          ROOT::Gdml
                          ROOT::FFTW
                          ROOT::Core
       )

target_compile_definitions(icarus_signal_processing_Filters INTERFACE EIGEN_FFTW_DEFAULT)

install_headers()
install_fhicl()
install_source()
