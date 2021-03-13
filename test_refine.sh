#!/bin/bash
set -xve
test_path=$1
radius=$2
test_case_name=$(basename $1)
outpath="testout/r$radius/$test_case_name/"
pb=$(realpath $test_path/*.pb)
vcf=$(realpath $test_path/*.vcf*)

renice -n 19 -p $$
mkdir -p $outpath
rm -f $outpath/*
{
build/usher --refine $radius -i $pb -v $vcf -o $outpath/refined.pb
} &> $outpath/log
/scratch/home/cheng/usher-master/build/usher -i $pb -v $vcf -o $outpath/not-refined.pb &> /dev/null
{
/scratch/home/cheng/usher-master/build/matUtils convert -i $outpath/refined.pb -v $outpath/refined.vcf
perl vcfNormalizeGt.pl $outpath/refined.vcf >$outpath/refined.normalized.vcf
}&
{
/scratch/home/cheng/usher-master/build/matUtils convert -i $outpath/not-refined.pb -v $outpath/not-refined.vcf
perl vcfNormalizeGt.pl $outpath/not-refined.vcf >$outpath/not-refined.normalized.vcf
}&
wait
{
diff -s $outpath/refined.normalized.vcf $outpath/not-refined.normalized.vcf
wc $outpath/refined.normalized.vcf
wc $outpath/not-refined.normalized.vcf
} &> $outpath/diff-out
