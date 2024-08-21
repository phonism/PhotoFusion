#!/bin/bash

function linux() {
    mkdir -p build
    cd build
    cmake \
        -DCMAKE_OSX_ARCHITECTURES="arm64;x86_64" \
        -DENABLE_LCMS=OFF \
        -DJPEG_INCLUDE_DIR=/usr/include \
        -DJPEG_LIBRARY=/usr/lib/x86_64-linux-gnu/libjpeg.so \
        -DBUILD_SHARED_LIBS=OFF \
        -DENABLE_JASPER=OFF \
        -DENABLE_EXAMPLES=OFF \
        -DUSE_OPENMP=OFF \
        -DUSE_FFTW=OFF \
        -DUSE_GCD=ON \
        -DCMAKE_VERBOSE_MAKEFILE=OFF \
        ../
    make -j 30
    cd ..

    mkdir -p output
    mkdir -p output/bin
    cp build/simple_raw2tiff output/bin/.
    cp build/simple_stack output/bin/.
}

function mac() {
    mkdir -p build
    cd build
    cmake \
        -DCMAKE_OSX_ARCHITECTURES="arm64;x86_64" \
        -DENABLE_LCMS=OFF \
        -DJPEG_INCLUDE_DIR=/Users/luqi03/workspace/github/libjpeg-turbo/install-arm64/include \
        -DJPEG_LIBRARY=/Users/luqi03/workspace/github/libjpeg-turbo/universal/libjpeg.a \
        -DBUILD_SHARED_LIBS=OFF \
        -DENABLE_JASPER=OFF \
        -DENABLE_EXAMPLES=OFF \
        -DUSE_OPENMP=OFF \
        -DUSE_FFTW=OFF \
        -DUSE_GCD=ON \
        -DCMAKE_VERBOSE_MAKEFILE=OFF \
        ../
    make -j 30
    cd ..

    mkdir -p output
    mkdir -p output/bin
    cp build/simple_raw2tiff output/bin/.
    cp build/simple_stack output/bin/.
}
mac
