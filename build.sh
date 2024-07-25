#!/bin/bash

mkdir -p build
cd build
cmake -DCMAKE_OSX_ARCHITECTURES="arm64;x86_64" ../
make -j 30
cd ..
