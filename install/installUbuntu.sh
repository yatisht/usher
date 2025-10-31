set -e
sudo -E apt update 
sudo -E apt-get install -yq --no-install-recommends build-essential \
 wget cmake  libboost-filesystem-dev libboost-program-options-dev libboost-iostreams-dev libboost-date-time-dev \
 libprotoc-dev libprotoc-dev protobuf-compiler \
 mafft rsync libtbb-dev openmpi-bin libopenmpi-dev automake libtool autoconf make nasm

# create build directory
startDir=$pwd
cd $(dirname "$0")
mkdir -p ../build
cd ../build

# download and build isa-l
wget https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz
tar -xvf v2.30.0.tar.gz
cd isa-l-2.30.0
./autogen.sh
./configure
make -j$(nproc)
sudo -E make install
cd ../..
echo ${pwd}

# build usher
cmake -S . -B build 
cmake --build build --parallel $(nproc)

# install faToVcf
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
chmod +x faToVcf

cd $startDir
