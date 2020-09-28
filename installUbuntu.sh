sudo apt update 
sudo apt-get --yes install build-essential
sudo apt-get --yes install wget
sudo apt-get --yes install cmake
sudo apt-get --yes install libboost-all-dev
sudo apt-get --yes install libomp-dev
sudo apt-get --yes install libprotoc-dev libprotoc-dev protobuf-compiler
sudo apt-get --yes install python python3-setuptools python3-pip 
sudo apt-get install -y mafft

wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_lin.tgz
tar -xvzf tbb2019_20191006oss_lin.tgz
mkdir -p build
cd build
cmake  -DTBB_DIR=${PWD}/../tbb2019_20191006oss  -DCMAKE_PREFIX_PATH=${PWD}/../tbb2019_20191006oss/cmake ..
make -j
cd ..
    
