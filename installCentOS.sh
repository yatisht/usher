sudo yum group install -y  "Development Tools"                                                                                                                                                                  
sudo yum install -y wget boost-devel rsync

# install cmake-3.18
wget https://github.com/Kitware/CMake/releases/download/v3.18.2/cmake-3.18.2.tar.gz
tar -xvzf cmake-3.18.2.tar.gz
cd cmake-3.18.2
./bootstrap --prefix=${PWD} --  -DCMAKE_USE_OPENSSL=OFF
make -j
make install
cd ..

# install mafft
wget https://mafft.cbrc.jp/alignment/software/mafft-7.471-gcc_fc6.x86_64.rpm
sudo rpm -Uvh mafft-7.471-gcc_fc6.x86_64.rpm
rm mafft-7.471-gcc_fc6.x86_64.rpm

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
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz 
tar -xvzf 2019_U9.tar.gz

# build programs
mkdir -p build
cd build
../cmake-3.18.2/bin/cmake  -DTBB_DIR=${PWD}/../oneTBB-2019_U9 -DTBB_ROOT=${PWD}/../oneTBB-2019_U9 -DCMAKE_PREFIX_PATH=${PWD}/../oneTBB-2019_U9/cmake  -DProtobuf_INCLUDE_DIRS=${PWD}/../protobuf-3.12.3/install/include/ -DProtobuf_LIBRARIES=${PWD}/../protobuf-3.12.3/cmake/build/libprotobuf.a -DProtobuf_PATH=${PWD}/../protobuf-3.12.3/cmake/build/lib64/cmake/protobuf ..
make -j
cd ..

# install faToVcf
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
chmod +x faToVcf
mv faToVcf build/
