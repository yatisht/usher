brew install cmake coreutils boost protobuf wget rsync openmpi libtool automake autoconf nasm isa-l

# create build directory
startDir=$pwd
cd $(dirname "$0")
cd ..

# TBB
#wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_mac.tgz
#tar -xvzf tbb2019_20191006oss_mac.tgz

# Build UShER
#cmake -DTBB_DIR=${PWD}/tbb2019_20191006oss -DCMAKE_PREFIX_PATH=${PWD}/tbb2019_20191006oss/cmake ..
#make -j2 VERBOSE=1
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17
cmake --build build --parallel 4

# install faToVcf
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/faToVcf .
chmod +x faToVcf

# install mafft
if ! command -v mafft &> /dev/null; then 
wget https://mafft.cbrc.jp/alignment/software/mafft-7.471-mac.zip
unzip mafft-7.471-mac.zip
cd mafft-mac/
mv mafft.bat /usr/local/bin/mafft; mv mafftdir /usr/local/bin/
cd ..
rm -rf mafft-mac
fi

cd $startDir
