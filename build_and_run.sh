#!/bin/bash
cd build
cmake ../cmake/local
cmake --build .
cd ..

./build/main
touch ./viewer/viewer.html

echo "finished"