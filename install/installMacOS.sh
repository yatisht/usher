brew install cmake coreutils boost protobuf wget rsync openmpi libtool automake autoconf nasm isa-l
brew upgrade protobuf

# create build directory
startDir=$pwd
cd $(dirname "$0")
mkdir -p ../build
cd ..

# build usher
cmake -S . -B build 
cmake --build build

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
