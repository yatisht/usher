sudo yum group install -y  "Development Tools"                                                                                                                                                                  
sudo yum install -y wget
sudo yum install -y boost-devel
sudo yum install -y python3


# install cmake-3.18
wget https://github.com/Kitware/CMake/releases/download/v3.18.2/cmake-3.18.2.tar.gz
tar -xvzf cmake-3.18.2.tar.gz
cd cmake-3.18.2
./bootstrap --prefix=${PWD} --  -DCMAKE_USE_OPENSSL=OFF
make -j
make install
cd ..

# install mafft
wget https://mafft.cbrc.jp/alignment/software/mafft-7.471-without-extensions-src.tgz
gunzip -cd mafft-7.471-without-extensions-src.tgz | tar xfv -
cd mafft-7.471-without-extensions/core/
make clean
make -j
sudo make install
cd ../../

# setup protobuf
wget https://github.com/protocolbuffers/protobuf/releases/download/v3.12.3/protobuf-cpp-3.12.3.tar.gz
tar -xvzf protobuf-cpp-3.12.3.tar.gz
cd protobuf-3.12.3
./configure --prefix=${PWD}/install
make -j2
make install
cd cmake; mkdir build; cd build;
../../../cmake-3.18.2/bin/cmake ..
make -j2
cd ../../../

# get TBB
wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_lin.tgz
tar -xvzf tbb2019_20191006oss_lin.tgz

# build programs
mkdir -p build
cd build
../cmake-3.18.2/bin/cmake  -DTBB_DIR=${PWD}/../tbb2019_20191006oss -DTBB_ROOT=${PWD}/../tbb2019_20191006oss -DCMAKE_PREFIX_PATH=${PWD}/../tbb2019_20191006oss/cmake  -DProtobuf_INCLUDE_DIRS=${PWD}/../protobuf-3.12.3/install/include/ -DProtobuf_LIBRARIES=${PWD}/../protobuf-3.12.3/cmake/build/libprotobuf.a -DProtobuf_PATH=${PWD}/../protobuf-3.12.3/cmake/build/lib64/cmake/protobuf ..
make -j
cd ..
