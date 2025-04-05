#!/bin/bash
set -e  # Exit immediately if any command fails

cmake -S . -B build
cmake --build build

./build/main
touch ./viewer/viewer.html

echo "finished"