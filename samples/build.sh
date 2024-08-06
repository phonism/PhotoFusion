#!/bin/bash

mkdir -p build
cd build
cmake \
    -DCMAKE_OSX_ARCHITECTURES="arm64;x86_64" \
    -DENABLE_LCMS=OFF \
    -DBUILD_SHARED_LIBS=OFF \
    -DENABLE_JASPER=OFF \
    -DENABLE_EXAMPLES=OFF \
    -DUSE_OPENMP=OFF \
    -DUSE_FFTW=OFF \
    -DUSE_GCD=ON \
    -DCMAKE_VERBOSE_MAKEFILE=OFF \
    -DJPEG_INCLUDE_DIR:PATH="" \
    -DJPEG_LIBRARY_DEBUG:FILEPATH="" \
    -DJPEG_LIBRARY_RELEASE:FILEPATH="" \
    ../
make -j 30
cd ..

mkdir -p output
mkdir -p output/bin
cp build/simple_raw2tiff output/bin/.
cp build/simple_stack output/bin/.
