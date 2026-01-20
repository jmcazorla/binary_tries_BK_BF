# if [ -d "build" ]; then
#     rm -rf build
# fi
# mkdir -p build/sdsl-install
mkdir build
cd build
cmake ..
make