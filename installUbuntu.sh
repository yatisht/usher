#sudo -E apt update 
#sudo -E apt-get --yes install build-essential \
# wget cmake  libboost-all-dev \
# libprotoc-dev libprotoc-dev protobuf-compiler \
# mafft rsync libtbb-dev

#download and install TBB
#wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz 
#tar -xvzf 2019_U9.tar.gz
mkdir -p build
cd build
cmake  -DTBB_DIR=${PWD}/../oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/../oneTBB-2019_U9/cmake ..
make -j4
cd ..
    
# install faToVcf
#rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
#chmod +x faToVcf
#mv faToVcf build/
