brew install cmake boost protobuf wget python@3

wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz 
tar -xvzf 2019_U9.tar.gz
mkdir -p build
cd build
cmake -DTBB_DIR=${PWD}/../oneTBB-2019_U9 -DCMAKE_PREFIX_PATH=${PWD}/../oneTBB-2019_U9/cmake ..
make -j
cd ..

# install mafft
if ! command -v mafft &> /dev/null; then 
wget https://mafft.cbrc.jp/alignment/software/mafft-7.471-mac.zip
unzip mafft-7.471-mac.zip
cd mafft-mac/
mv mafft.bat /usr/local/bin/mafft; mv mafftdir /usr/local/bin/
cd ..
rm -rf mafft-mac
fi
