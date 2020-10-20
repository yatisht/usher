sudo apt update 
sudo apt-get --yes install build-essential
sudo apt-get --yes install wget
sudo apt-get --yes install cmake
sudo apt-get --yes install libboost-all-dev
sudo apt-get --yes install libprotoc-dev libprotoc-dev protobuf-compiler
sudo apt-get --yes install python python3-setuptools python3-pip 
sudo apt-get install -y mafft

wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz 
tar -xvzf 2019_U9.tar.gz
mkdir -p build
cd build
cmake  -DTBB_DIR=${PWD}/../oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/../oneTBB-2019_U9/cmake ..
make -j
cd ..
    
# install faToVcf
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToVcf
chmod 777 ./faToVcf
mv ./faToVcf ./scripts/faToVcf
