#!/bin/bash -ex
#
# Get the raw sequences for all the descendents

# Raw sequence file was copied from GCP Storage Bucket to local directory, 
# passed to this script through first argument.
all_sequences_fasta=$1

# SARS-CoV-2 reference genome copied from GCP Storage Bucket to local directory, 
# passed to this script through second argument.
reference=$2  
startDir=$PWD
cores=`grep -c ^processor /proc/cpuinfo`

# Check correct number of args passed
if [ "$#" -ne 2 ]; then
    echo "ERROR: Incorrect number of arguments passed."
fi

mkdir -p filtering/fastas
cp $reference filtering/fastas/reference.fa
cp $all_sequences_fasta filtering/fastas/extractedSeqs.fa

mkdir -p filtering/fastas/OrderedRecombs
mkdir -p filtering/fastas/AlignedRecombs
python3 filtering/analyzerecomb.py -a

# Align raw sequences using mafft
cd filtering/fastas/OrderedRecombs
ls . |  parallel -j $cores "mafft --auto {} > ../AlignedRecombs/{} "
cd $startDir

# Generate report
cp filtering/empty_report.txt filtering/data/report.txt
ls filtering/fastas/AlignedRecombs | parallel -j $cores python3 filtering/checkmutant.py {.} -r 

awk '$19 == "False"' filtering/data/report.txt | awk '$14 == "False"' | awk '$11 == "False"' > filtering/data/final_report.txt
echo "DONE"
