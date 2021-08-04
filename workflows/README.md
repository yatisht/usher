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

RIPPLES: detect recombinants in the ancestry of the user-supplied samples

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="ripples"

## Further Reading:

More information about each of these utilities can be found here: https://usher-wiki.readthedocs.io/en/latest/
