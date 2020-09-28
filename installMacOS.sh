brew install cmake boost libomp protobuf wget python@3

wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_mac.tgz
tar -xvzf tbb2019_20191006oss_mac.tgz
mkdir -p build
cd build
cmake -DTBB_DIR=${PWD}/../tbb2019_20191006oss -DCMAKE_PREFIX_PATH=${PWD}/../tbb2019_20191006oss/cmake -DOpenMP_CXX_FLAGS="-Xclang -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_C_FLAGS="-Xclang -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_CXX_LIB_NAMES=libomp  -DOpenMP_C_LIB_NAMES=libomp  -DOpenMP_libomp_LIBRARY=/usr/local/opt/libomp/lib/libomp.dylib -DCMAKE_SHARED_LINKER_FLAGS="-L/usr/local/opt/libomp/lib -lomp -Wl,-rpath,/usr/local/opt/libomp/lib" ..
make -j
cd ..

# install mafft
wget https://mafft.cbrc.jp/alignment/software/mafft-7.471-mac.zip
unzip mafft-7.471-mac.zip
cd mafft-mac/
mv mafft.bat mafft
export PATH=$PATH:$(pwd)
