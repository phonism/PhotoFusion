#!/bin/bash

cd build
cmake -DCMAKE_OSX_ARCHITECTURES=x86_64 ../
make -j 30
cd ..

mkdir -p output
mkdir -p output/bin
cp build/simple_raw2tiff output/bin/.
