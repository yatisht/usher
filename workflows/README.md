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

outbreak: run extract on the dataset that includes user provided samples to identify close related sequences

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="outbreak"

augur: runs the augur pipeline to build a clocked tree that includes the user samples for visualization

    snakemake --use-conda --cores [num threads] --config FASTA="[user_fa]" RUNTYPE="augur"
    
Note that adding "-d [run_dir]" to the command line above will generate all output files in the specified directory. To do this, you must provide the full path to the fasta file or place the fasta file into the specified run directory. 

## Quick Tutorial:

To run a file through any of the workflows, with the exception of _introduce_, replace _user_fa_ in the command line with "_data/example.fasta_". To run a file through the _introduce_ workflow, use "_data/introduce_example.fasta_".

### Augur workflow

The _augur_ workflow produces a json file with published sequences provided by the user that belong to a certain outbreak. The user should alter the _config.json_ file within the config folder to describe their dataset, as the current file is set as an example of the workflow.

For further help with the _Augur_ workflow please contact [Adriano Schneider](mailto:adeberna@ucsc.edu).

## Further Reading:

More information about each of these utilities can be found here: https://usher-wiki.readthedocs.io/en/latest/
