#!/bin/bash

# Navigate to the build directory
cd build

echo "Clean the build directory (remove all files in the directory)"
rm -rf *

echo "Run cmake and build the project"
cmake ../cmake/local
cmake --build .

echo "Finished cmake build"
