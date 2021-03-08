sudo -E apt update 
sudo -E apt-get --yes install build-essential \
 wget cmake  libboost-all-dev \
 libprotoc-dev libprotoc-dev protobuf-compiler \
 python python3-setuptools python3-pip mafft rsync
 
#install python library biopython
pip3 install biopython

#download and install TBB
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz 
tar -xvzf 2019_U9.tar.gz
mkdir -p build
cd build
cmake  -DTBB_DIR=${PWD}/../oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/../oneTBB-2019_U9/cmake ..
make -j
cd ..
    
# install faToVcf
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
chmod +x faToVcf
