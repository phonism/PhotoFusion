#!/bin/bash

set -eu

# Values set in env-${PLATFORM}.sh will override these settings
# Modify env-${PLATFORM}.sh to change these settings, not here
#CXX_COMPILER=/usr/local/opt/llvm/bin/clang++
CXX_COMPILER=clang++
BUILD_CONFIG=Release
BASE_DIR=./
PLATFORM=iOS
CMAKE_SYSTEM_NAME=iOS
CMAKE_OSX_SYSROOT=iphoneos
BINARY_PACKAGE_BUILD=ON
ARCHS=arm64

CONFIG_FILE="${BASE_DIR}/env-${PLATFORM}.sh"

if [ -f "${CONFIG_FILE}" ]; then
  . "${CONFIG_FILE}"
fi

SOURCE_DIR="${BASE_DIR}/."

BUILD_DIR="${BASE_DIR}/build-${PLATFORM}"

config_and_build() {
  rm -rf "${BUILD_DIR}"

  cmake \
  -DCMAKE_SYSTEM_NAME=$CMAKE_SYSTEM_NAME \
  -DCMAKE_OSX_SYSROOT=$CMAKE_OSX_SYSROOT \
  -DCMAKE_OSX_ARCHITECTURES="${ARCHS}" \
    -DJPEG_INCLUDE_DIR=/Users/luqi03/workspace/github/libjpeg-turbo/install-arm64/include \
    -DJPEG_LIBRARY=/Users/luqi03/workspace/github/libjpeg-turbo/universal/libjpeg.a \
  -DBUILD_SHARED_LIBS=OFF \
  -DENABLE_LCMS=OFF \
  -DENABLE_JASPER=OFF \
  -DUSE_OPENMP=OFF \
  -DUSE_JPEG=OFF \
  -DUSE_FFTW=OFF \
  -DUSE_GCD=ON \
  -DENABLE_EXAMPLES=OFF \
  -DCMAKE_VERBOSE_MAKEFILE=OFF \
  -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE \
  -DCMAKE_BUILD_TYPE:STRING=$BUILD_CONFIG \
  -DCMAKE_CXX_COMPILER:FILEPATH="${CXX_COMPILER}" \
  -DCMAKE_CXX_STANDARD_LIBRARIES:STRING="-stdlib=libc++" \
  -S"${SOURCE_DIR}" \
  -B"${BUILD_DIR}" \
  -G Ninja

  cmake \
  --build "${BUILD_DIR}" \
  --config $BUILD_CONFIG \
  --verbose \
  --target PhotoFusionAPI --
}
config_and_build
