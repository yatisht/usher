# This directory contains simple snakemake workflows for running UShER, RIPPLES, and matUtils

## To use these workflows, include in the current working directory:

1. a fasta file with SARS-CoV-2 genome sequences: [user_fa]
2. the Snakefile
3. the conda environment file, usher.yaml

## Users can run each workflow as:

UShER: add samples to the latest public MAT

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="usher"

matUtils: extract subtrees in auspice.us compatible json format using matUtils

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="matUtils"
    
taxodium: view the "big tree" in the taxodium (taxodium.org) tree viewing platform

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="taxodium" 
    
translate: output lineage-aware translation for each amino acid substitution

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="translate" 

RIPPLES: detect recombinants in the ancestry of the user-supplied samples

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="ripples"
    
introduce: search for unique introductions within the user-supplied samples

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="introduce"

systematic: search for possible systematic errors in your added samples by outputing a list of sites whose parsimony score increased

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="systematic"

Note that adding "-d [run_dir]" to the command line above will generate all output files in the specified directory. To do this, you must provide the full path to the fasta file or place the fasta file into the specified run directory. 

## Further Reading:

More information about each of these utilities can be found here: https://usher-wiki.readthedocs.io/en/latest/
