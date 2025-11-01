#brew install cmake coreutils boost protobuf wget rsync openmpi libtool automake autoconf nasm isa-l

# Check if conda is available
if ! command -v conda >/dev/null 2>&1; then
    echo "WARNING: Installing conda since it is not installed or not found in your PATH." >&2
    curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o miniconda.sh  
    bash miniconda.sh -b -p "$HOME/miniconda" 
    rm miniconda.sh 
    export PATH="$HOME/miniconda/bin:$PATH" 
    source "$HOME/miniconda/etc/profile.d/conda.sh"
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
if ! command -v mafft >/dev/null 2>&1; then
    wget https://mafft.cbrc.jp/alignment/software/mafft-7.471-mac.zip
    unzip mafft-7.471-mac.zip
    cd mafft-mac/
    mv mafft.bat /usr/local/bin/mafft; mv mafftdir /usr/local/bin/
    cd ..
    rm -rf mafft-mac
fi

cd $startDir
