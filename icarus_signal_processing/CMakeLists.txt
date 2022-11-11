add_subdirectory(Filters)
add_subdirectory(Detection)


FILE(GLOB src_files *.cxx)

cet_make_library( SOURCE ${src_files}
          HEADERS_TARGET
          LIBRARIES
                          icarus_signal_processing_Detection
                          icarus_signal_processing_Filters
                          PkgConfig::FFTW3F
                          PkgConfig::FFTW3
                          PkgConfig::FFTW3L
                          PkgConfig::FFTW3Q
                          ROOT::Geom
                          ROOT::XMLIO
                          ROOT::Gdml
                          ROOT::FFTW
                          ROOT::Core
       )

target_compile_definitions(icarus_signal_processing INTERFACE EIGEN_FFTW_DEFAULT)

install_headers()
install_fhicl()
install_source()
