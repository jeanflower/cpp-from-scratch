#!/bin/bash
set -e  # Exit immediately if any command fails

cmake -S . -B build
cmake --build build

cd build
ctest -V -O ../output/test_output.log
