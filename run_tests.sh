#!/bin/bash

cd cmake_build
ctest --verbose -D ExperimentalMemCheck .
cd ..
