cd build

echo "Clean the build directory (remove all files in the directory)"
rm -rf *

cd ..

clang++ -g -o build/main main.cpp

# will generate binary file called main
# will generate debugging symbols as main.dSYM

echo "Finished clang build"
