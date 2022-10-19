add_subdirectory(Filters)
add_subdirectory(Detection)

add_definitions(-DEIGEN_FFTW_DEFAULT)
include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})

cet_make_library( SOURCE 
                          Denoising.cxx
                          ROIFinder2D.cxx
                          WaveformTools.cxx
          HEADERS_TARGET
          LIBRARIES 
                          icarus_signal_processing_Detection
                          icarus_signal_processing_Filters
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
