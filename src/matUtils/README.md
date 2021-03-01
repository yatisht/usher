# matUtils
matUtils is a toolkit for working with mutation annotated tree protobuf files, targeted towards filtering and extraction of desired information. It is divided into a number of subcommands. Each subcommand is called individually and has its own parameters and options; for example "matUtils annotate --help" yields the help message for the annotate subcommand of the matUtils toolkit.

### Common Options

**-i**: Input mutation-annotated tree file (REQUIRED)

**-T**: Number of threads to use for multithreaded procedures. Default is all available

**-h**: Print help messages

## Annotate

The annotate command takes a MAT protobuf file and a two-column tab-separated text file indicating sample names and lineage assignments. The software will automatically identify the best clade root for that lineage and save the assignment to each sample indicated, returning a new protobuf with these values stored.

### Specific Options

**-o**: Output mutation-annotated tree file (REQUIRED)

**-l**: Provide a path to a file assigning lineages to samples to locate and annotate clade root nodes (REQUIRED)

**-f**: Set the minimum allele frequency in samples to find the best clade root node (default = 0.9)

## Convert 

The convert subcommand yields non-protobuf format files representing the tree. This can be a newick format file and/or a VCF representing each sample's variations from the reference. 

### Specific Options

**-v**: Output VCF file 

**-t**: Output Newick tree file

**-n**: Do not include sample genotype columns in VCF output. Used only with -v

## Mask

The mask subcommand restricts indicated sample names and returns a masked protobuf.

### Specific Options

**-o**: Output mutation-annotated tree file (REQUIRED)

**-s**: Use to mask specific samples from the tree. 

## Describe

Describe fetches the path of mutations associated with each indicated sample's placement on the tree and prints the results.

### Specific Options
 
**-m**: File containing sample names for which mutation paths should be displayed (REQUIRED)

## Filter

Filter removes an indicated set of samples from the input MAT, returning a smaller MAT.

### Specific Options

**-o**: Output mutation annotated tree [REQUIRED]

**-p**: File containing names of samples to be pruned from the input MAT.

**-P**: File containing names of samples to retain in the input MAT, pruning all others.

## Uncertainty

The uncertainty command calculates placement quality metrics for a set of samples or the whole tree and returns these values in the indicated format.

### Specific Options

**-s**: File containing names of samples to calculate metrics for.

**-g**: Calculate and print the total tree parsimony score. 

**-e**: Name for an output Nextstrain Auspice-compatible .tsv file of the number of equally parsimonious placement values for each indicated sample- smaller is better, best is 1 (requires -s).

**-n**: Name for an output Nextstrain Auspice-compatible .tsv file of neighborhood size scores for the equally parsimonious placements for each indicated sample- smaller is better, best is 0 (requires -s).

## Summary

This command gets basic statistics about the input MAT.

### Specific Options

**-s**: Write a tsv listing all samples in the tree and their parsimony scores (terminal branch length).

**-c**: Write a tsv listing all clades in the tree and their occurrence over nodes in the tree.

**-m**: Write a tsv listing all mutations in the tree and their occurrence count.
