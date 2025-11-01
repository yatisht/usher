#brew install cmake coreutils boost protobuf wget rsync openmpi libtool automake autoconf nasm isa-l

# Check if conda is available
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda is not installed or not found in your PATH." >&2
    echo "Please install Miniconda or Anaconda before proceeding." >&2
    exit 1
fi

# create build directory
startDir=$pwd
cd $(dirname "$0")
conda env create -f environment.yml -n usher -y
conda activate usher
cd ..

# build usher
cmake -S . -B build 
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
