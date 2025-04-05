#!/bin/bash

# Navigate to the build directory
cd build

echo "Clean the build directory (remove all files in the directory)"
rm -rf *

cd ..

echo "Run cmake and build the project"

cmake -S . -B build
cmake --build build

echo "Finished cmake build"
