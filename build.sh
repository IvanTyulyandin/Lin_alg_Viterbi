#!/bin/bash

mkdir -p cmake_build
cd cmake_build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_CXX_COMPILER=clang++ \
    ..
make -j 5
cd ..
