# matUtils
matUtils is a toolkit for working with mutation annotated tree protobuf files, targeted towards filtering and extraction of desired information. It is divided into a number of subcommands. Each subcommand is called individually and has its own parameters and options; for example "matUtils annotate --help" yields the help message for the annotate subcommand of the matUtils toolkit.

### Common Options

**-i**: Input mutation-annotated tree file (REQUIRED)

**-T**: Number of threads to use for multithreaded procedures. Default is all available

**-h**: Print help messages

## Extract

This command selects subtrees from the MAT based on clade identification, sample choice files, or placement quality metrics, and optionally converts to other file formats, including vcf and newick.

**-s**: Select samples to extract by explicitly naming them, one per line in a text file.

**-c**: Select samples by membership in at least one of the indicated clade(s), comma delimited.

**-m**: Select samples by whether they contain any of the indicated mutation(s), comma delimited.

**-e**: Select samples to extract by whether they have less than the indicated number of equally parsimonious placements. Note: adds significantly to runtime with large sample selections.

**-a**: Select samples to extract by whether their parsimony score (terminal branch length) is less than the indicated value.

**-r**: Automatically select two samples per unique clade identifier to extract.

**-p**: Prune the samples indicated by other arguments instead of extracting them to create the subtree.

**-d**: Set the directory to save output files. Directory must already exist. Default is current directory

**-S**: Write the path of mutations defining each selected sample to a file.

**-C**: Write the path of mutations defining each clade in the subtree after sample extraction to a file.

**-A**: Write mutations assigned to each node in the selected tree to the target file.

**-v**: Convert the subtree to a VCF and write it to the indicated file.

**-n**: Do not include sample genotype columns in the VCF output. Used only with -v

**-o**: Create a new MAT .pb file representing the extracted subtree. 

**-O**: Collapse the subtree before saving.

**-t**: Write a newick string representing the extracted subtree to the indicated file.

## Summary

This command gets basic statistics about the input MAT.

### Specific Options

**-s**: Write a tsv listing all samples in the tree and their parsimony scores (terminal branch length).

**-c**: Write a tsv listing all clades in the tree and their occurrence over nodes in the tree.

**-m**: Write a tsv listing all mutations in the tree and their occurrence count.

## Annotate

The annotate command takes a MAT protobuf file and a two-column tab-separated text file indicating sample names and lineage assignments. The software will automatically identify the best clade root for that lineage and save the assignment to each sample indicated, returning a new protobuf with these values stored. Optionally, it can take a tsv directly mapping clades to internal node identifiers (with -C) and assign labels accordingly without inference.

### Specific Options

**-o**: Output mutation-annotated tree file (REQUIRED)

**-l**: Provide a path to a file assigning lineages to samples to locate and annotate clade root nodes (REQUIRED)

**-f**: Set the minimum allele frequency in samples to find the best clade root node (default = 0.9)

**-c**: Path to a file containing clade asssignments of samples. An algorithm automatically locates and annotates clade root nodes.

**-C**: Path to a tsv file mapping clades to their respective internal node identifiers. Used to directly assign labels without inference.

**-s**: Minimum fraction of the clade samples that should be desecendants of the assigned clade root (default = 0.6)

**-l**: Clear current clade annotations before applying new ones.

## Uncertainty

The uncertainty command calculates placement quality metrics for a set of samples or the whole tree and returns these values in the indicated format.

### Specific Options

**-s**: File containing names of samples to calculate metrics for.

**-g**: Calculate and print the total tree parsimony score. 

**-e**: Name for an output Nextstrain Auspice-compatible .tsv file of the number of equally parsimonious placement values for each indicated sample- smaller is better, best is 1 (requires -s).

**-n**: Name for an output Nextstrain Auspice-compatible .tsv file of neighborhood size scores for the equally parsimonious placements for each indicated sample- smaller is better, best is 0 (requires -s).

## Mask [EXPERIMENTAL]

The mask subcommand restricts indicated sample names and returns a masked protobuf.

### Specific Options

**-o**: Output mutation-annotated tree file (REQUIRED)

**-s**: Use to mask specific samples from the tree. 
