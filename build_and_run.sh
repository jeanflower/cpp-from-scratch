#!/bin/bash
set -e  # Exit immediately if any command fails

cmake -S . -B build
cmake -DCMAKE_BUILD_TYPE=Release -B build
#cmake -DCMAKE_BUILD_TYPE=Debug -B build
cmake --build build --target main

./build/main
touch ./viewer/viewer.html

echo "finished"