# - Try to find the FFTW library
# Once done this will define
#  FFTW_FOUND - system has FFTW
#  FFTW_INCLUDE_DIRS - the FFTW include directory
#  FFTW_LIBRARIES - link these to use FFTW

find_path(FFTW_INCLUDE_DIR fftw3.h
  HINTS
    $ENV{FFTW_DIR}/include
    /usr/local/include
    /usr/include
)

find_library(FFTW_LIBRARY NAMES fftw3
  HINTS
    $ENV{FFTW_DIR}/lib
    /usr/local/lib
    /usr/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDE_DIR FFTW_LIBRARY)

if(FFTW_FOUND)
  set(FFTW_LIBRARIES ${FFTW_LIBRARY})
  set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
endif()