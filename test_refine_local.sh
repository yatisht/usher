#!/bin/bash
mkdir -p testout
rm -f testout/*
build/usher --refine $3 -i $1 -v $2 -o testout/refined.pb
../usher-master/build/usher -i $1 -v $2 -o testout/not-refined.pb
mkfifo testout/refined.vcf
mkfifo testout/not-refined.vcf
../usher-master/build/matUtils convert -i testout/refined.pb -v testout/refined.vcf&
perl vcfNormalizeGt.pl testout/refined.vcf >testout/refined.normalized.vcf
../usher-master/build/matUtils convert -i testout/not-refined.pb -v testout/not-refined.vcf&
perl vcfNormalizeGt.pl testout/not-refined.vcf >testout/not-refined.normalized.vcf
diff -s testout/refined.normalized.vcf testout/not-refined.normalized.vcf
wc testout/refined.normalized.vcf
wc testout/not-refined.normalized.vcf
