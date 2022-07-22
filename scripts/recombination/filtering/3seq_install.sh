#!/bin/bash -ex
#
# Install 3seq 
startDir=$PWD
git clone https://gitlab.com/mrkylesmith/3seq.git
cd 3seq
make
./3seq -g my3seqTable700 700
cd $startDir
