# Ultrafast Sample Placement on Existing Trees (UShER)

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/yatisht/usher/blob/master/LICENSE

[![License][license-badge]][license-link]
[![Build Status](https://github.com/yatisht/usher/workflows/build/badge.svg)](https://github.com/yatisht/usher/actions)

<img src="/images/usher_logo.png" width="600">

**NEW: We will now be sharing and updating UShER's pre-processed mutation-annotated tree object for public SARS-CoV-2 sequences here: https://hgwdev.gi.ucsc.edu/~angie/UShER_SARS-CoV-2/**

UShER is a program that rapidly places new samples onto an existing phylogeny using maximum parsimony. It is particularly helpful in understanding the relationships of newly sequenced SARS-CoV-2 genomes with each other and with previously sequenced genomes in a global phylogeny. This has emerged as an important challenge during this pandemic for enabling *genomic contact tracing* since the viral phylogeny is already very large (>100K sequences, see https://github.com/roblanf/sarscov2phylo/releases) and is expected to grow by many fold in the coming months. 

UShER is much faster than existing tools with similar functionality and has now also been integrated in the [UCSC SARS-CoV-2 Genome Browser](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace), which does not require UShER installation and usage know-how as described below for SARS-CoV-2 applications. Please follow the steps below if you wish to use UShER in a standalone fashion.

## Contents
* [Installing UShER](#installing-usher)
    - [Using Docker](#using-docker)
        - [Local build and install](#local-build-and-install)
        - [Using pre-built Docker image](#using-pre-built-docker-image)
    - [Using Conda](#using-conda)
    - [Using installation scripts](#using-installation-scripts)
* [How UShER works?](#how-usher-works)
* [Using UShER](#using-usher)
  - [Displaying help message](#displaying-help-message)
  - [Pre-processing global phylogeny](#pre-processing-global-phylogeny)
  - [Placing new samples](#placing-new-samples)
  - [Uncertainty in placing new samples](#uncertainty-in-placing-new-samples)
      * [Branch Parsimony Score](#branch-parsimony-score)
      * [Multiple parsimony-optimal placements](#multiple-parsimony-optimal-placements)
      * [Updating multiple input trees](#updating-multiple-input-trees)
* [Fasta2UShER](#fasta2usher)
* [matUtils](#matUtils)
* [Acknowledgement](#acknowledgement)
* [Reference](#reference)

## Installing UShER

First, clone the UShER repository using the commands below.
```
git clone https://github.com/yatisht/usher.git
cd usher
```

Next, you may install UShER either using [Docker](https://www.docker.com/) (**recommended**) or using one of the installation scripts that we provide for **MacOS (10.14 and above)**, **Ubuntu (18.04 and above)** and **CentOS (7 and above)**. See instructions below.

### Using Docker

Perhaps the simplest method to install this package irrespective of the platform (Windows, Linux or MacOS) is using [Docker](https://www.docker.com/) (see Docker installation instructions [here](https://docs.docker.com/get-docker/)). 

#### Local build and install

The following command to build and generate your own UShER image locally.

```
docker build --no-cache -t usher .
```
Once the image is generated, you may start a new bash in a Ubuntu 18.04 container using the Docker command below, following which the remaining commands on this page should work smoothly.
```
docker run -t -i usher /bin/bash
```

#### Using pre-built Docker image

Alternatively, you can also use a pre-built Docker image shared by us using the following commands.

```
docker pull yatisht/usher:latest
docker run -t -i yatisht/usher:latest /bin/bash
```

#### Using Conda
Installing UShER:
```
conda env create -f environment.yml
conda activate usher
git clone https://github.com/oneapi-src/oneTBB
cd oneTBB
git checkout cc2c04e2f5363fb8b34c10718ce406814810d1e6
cd ..
mkdir build
cd build
cmake  -DTBB_DIR=${PWD}/../oneTBB  -DCMAKE_PREFIX_PATH=${PWD}/../oneTBB/cmake ..
make -j
cd ..
``` 
[Fasta2UShER](#fasta2usher) requires the faToVcf utility that can be obtained as follows.

* MacOS
```
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/faToVcf .
chmod +x faToVcf
mv faToVcf scripts/
```
* Linux
```
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
chmod +x faToVcf
mv faToVcf scripts/
```

### Using installation scripts

If you don't have Docker/Conda or prefer not to use it for some reason, you may use one of our installation scripts below to install UShER.

* For MacOS 10.14 and above: 
```
./installMacOS.sh  
```
* For Ubuntu 18.04 and above (requires sudo privileges):
```
./installUbuntu.sh  
```
* For CentOS 7 and above (requires sudo privileges): 
```
./installCentOS.sh  
```

## How UShER works?

Given existing samples, whose genotypes and phylogenetic tree is known, and the genotypes of new samples, UShER aims to incorporate new samples into the phylogenetic tree while preserving the topology of existing samples and maximizing parsimony. UShER’s algorithm consists of two phases: (i) the pre-processing phase and (ii) the placement phase.

In the **pre-processing phase**, UShER accepts the phylogenetic tree of existing samples in a Newick format and their genotypes, specified as a set of single-nucleotide variants with respect to a reference sequence (UShER currently ignores indels), in a VCF format. For each site in the VCF, UShER uses [Fitch-Sankoff algorithm](https://evolution.gs.washington.edu/gs541/2010/lecture1.pdf) to find the most parsimonious nucleotide assignment for every node of the tree (UShER automatically labels internal tree nodes). When a sample contains **ambiguous genotypes**, multiple nucleotides may be most parsimonious at a node. To resolve these, UShER assigns it any one of the most parsimonious nucleotides with preference, when possible, given to the reference base. UShER also allows the VCF to specify ambiguous bases in samples using [IUPAC format](https://www.bioinformatics.org/sms/iupac.html), which are also resolved to a unique base using the above strategy. When a node is found to carry a mutation, i.e. the base assigned to the node differs from its parent, the mutation gets added to a list of mutations corresponding to that node. Finally, UShER uses [protocol buffers](https://developers.google.com/protocol-buffers) to store in a file, the Newick string corresponding to the input tree and a list of lists of node mutation, which we refer to as **mutation-annotated tree object**, as shown in the figure below.

![pre-processing](/images/pre-processing.png)

The mutation-annotated tree object carries sufficient information to derive parsimony-resolved genotypes for any tip of the tree using the sequence of mutations from the root to that tip. For example, in the above figure, S5 can be inferred to contain variants G1149U, C7869U, G3179A and A2869G with respect to the reference sequence. Compared to other tools that use full multiple-sequence alignment (MSA) to guide the placement, UShER's mutation-annotated tree object is compact and is what helps make it **fast**.

In the **placement phase**, UShER loads the pre-processed mutation-annotated tree object and the genotypes of new samples in a VCF format and **sequentially** adds the new samples to the tree. For each new sample, UShER computes the additional parsimony score required for placing it at every node in the current tree while considering the full path of mutations from the root of the tree to that node. Next, UShER places the new sample at the node that results in the smallest additional parsimony score. When multiple node placements are equally parsimonious, UShER picks the node with a greater number of descendant leaves for placement. If the choice is between a parent and its child node, the parent node would always be selected by this rule. However, a more accurate placement should reflect the number of leaves uniquely attributable to the child versus parent node. Therefore, in these cases, UShER picks the parent node if the number of descendant leaves of the parent that are not shared with the child node exceed the number of descendant leaves of the child. The figure below shows a new sample, S7, containing variants G1149U and C9977A being added to the previous mutation-annotated tree object in a parsimony-optimal fashion (with a parsimony score of 1 for the mutation C9977A). UShER also automatically imputes and reports **ambiguous genotypes** for the newly added samples and ignores **missing bases**, such as 'N' or '.' (i.e. missing bases never contribute to the parsimony score).

![placement](/images/placement.png)

At the end of the placement phase, UShER allows the user to create another protocol-buffer (protobuf) file containing the mutation-annotated tree object for the newly generated tree including added samples as also shown in the example figure above. This allows for another round of placements to be carried out over and above the newly added samples. 

## Using UShER

### Displaying help message

To familiarize with the different command-line options of UShER, it would be useful to view its help message using the command below:

```
./build/usher --help
```

### Pre-processing global phylogeny

The following example command pre-processes the existing phylogeny (`global_phylo.nh`) and using the genotypes (`global_samples.vcf`) and generates the mutation-annotated tree object that gets stored in a protobuf file (`global_assignments.pb`). Note that UShER would automatically place onto the input global phylogeny any samples in the VCF (to convert a fasta sequence to VCF, consider using [Fasta2USHER](#fasta2usher)) that are missing in the input global phylogeny using its parsimony-optimal placement algorithm. This final tree is written to a file named `final-tree.nh` in the folder specified by `--outdir` or `-d` option (if not specified, default uses current directory). 
```
./build/usher -t test/global_phylo.nh -v test/global_samples.vcf -o global_assignments.pb -d output/
```
By default, UShER uses **all available threads** but the user can also specify the number of threads using the `--threads` or `-T` command-line parameter.

UShER also allows an option during the pre-processing phase to collapse nodes (i.e. delete the node after moving its child nodes to its parent node) that are not inferred to contain a mutation through the Fitch-Sankoff algorithm as well as to condense nodes that contain identical sequences into a single representative node. This is the **recommended usage** for UShER as it not only helps in significantly reducing the search space for the placement phase but also helps reduce ambiguities in the placement step and can be done by setting the `--collapse-tree` or `-c` parameter. The collapsed input tree is stored as `condensed-tree.nh` in the output directory. 
```
./build/usher -t test/global_phylo.nh -v test/global_samples.vcf -o global_assignments.pb -c -d output/
```

Note the the above command would condense identical sequences, namely S2, S3 and S4, in the example figure above into a single condensed new node (named something like *node_1_condensed_3_leaves*). If you wish to display the collapsed tree without condensing the nodes, also set the `--write-uncondensed-final-tree` or `-u` option, for example, as follows:
```
./build/usher -t test/global_phylo.nh -v test/global_samples.vcf -o global_assignments.pb -c -u -d output/
```
The above commands saves the collapsed but uncondensed tree as `uncondensed-final-tree.nh` in the output directory. 

### Placing new samples

Once the pre-processing is complete and a mutation-annotated tree object is generate (e.g. `global_assignments.pb`), UShER can place new sequences whose variants are called in a VCF file (e.g. `new_samples.vcf`) to existing tree as follows:

```
./build/usher -i global_assignments.pb -v test/new_samples.vcf -u -d output/
```
Again, by default, UShER uses **all available threads** but the user can also specify the number of threads using the *--threads* command-line parameter.

The above command not only places each new sample sequentially, but also reports the parsimony score and the number of parsimony-optimal placements found for each added sample. UShER displays warning messages if several (>=4) possibilities of parsimony-optimal placements are found for a sample. This can happen due to several factors, including (i) missing data in new samples, (ii) presence of ambiguous genotypes in new samples and (iii) structure and mutations in the global phylogeny itself, including presence of multiple back-mutations. 

In addition to the global phylogeny, one often needs to contextualize the newly added sequences using subtrees of closest *N* neighbouring sequences, where *N* is small. UShER allows this functionality using `--write-subtrees-size` or `-k` option, which can be set to an arbitrary *N*, such as 20 in the example below:

```
./build/usher -i global_assignments.pb -v test/new_samples.vcf -u -k 20 -d output/
```
The above command writes subtrees to files names `subtree-<subtree-number>.nh`. It also write a text file for each subtree (named `subtree-<subtree-number>-mutations.txt` showing mutations at each internal node of the subtree. If the subtrees contain condensed nodes, it writes the expanded leaves for those nodes to text files named `subtree-<subtree-number>-expanded.txt`. 

Finally, the new mutation-annotated tree object can be stored again using `--save-mutation-annotated-tree` or `-o` option (overwriting the loaded protobuf file is allowed).

```
./build/usher -i global_assignments.pb -v test/new_samples.vcf -u -o new_global_assignments.pb -d output/
```


### Uncertainty in placing new samples

#### Branch Parsimony Score

UShER also allows quantifying the uncertainty in placing new samples by reporting the parsimony scores of adding new samples to all possible nodes in the tree **without** actually modifying the tree (this is because the tree structure, as well as number of possible optimal placements could change with each new sequential placement). In particular, this can help the user explore which nodes of the tree result in a small and optimal or near-optimal parsimony score. This can be done by setting the `--write-parsimony-scores-per-node` or `-p` option, for example, as follows:
```
./build/usher -i global_assignments.pb -v test/new_samples.vcf -p -d output/
```
The above command writes a file `parsimony-scores.tsv` containing branch parsimony scores to the output directory. Note that because the above command does not perform the sequential placement on the tree, the number of parsimony-optimal placements reported for the second and later samples could differ from those reported with actual placements.

The figure below shows how branch parsimony score could be useful for uncertainty analysis. The figure shows color-coded parsimony score of placing a new sample at different branches of the tree with black arrow pointing to the branch where the placement is optimal. As can be seen from the color codes, the parsimony scores are low (implying good alternative placement) for several neighboring branches of the optimal branch. 

<img src="/images/branch-parsimony-score.png" width="400">

#### Multiple parsimony-optimal placements

To further aid the user to quantify phylogenetic uncertainty in placement, UShER has an ability to enumerate all possible topologies resulting from equally parsimonious sample placements. UShER does this by maintaining a list of mutation-annotated trees (starting with a single mutation-annotated tree corresponding to the input tree of existing samples) and sequentially adds new samples to each tree in the list while increasing the size of the list as needed to accommodate multiple equally parsimonious placements for a new sample. This feature is available using the `--multiple-placements` or `-M` option in which the user specifies the maximum number of topologies that UShER should maintain before it reverts back to using the default tie-breaking strategy for multiple parsimony-optimal placements in order to keep the runtime and memory usage of UShER reasonable. 

```
./build/usher -i global_assignments.pb -v <USER_PROVIDED_VCF> -M -d output/
```

Note that if the number of equally parsimonious placements for the initial samples is large, the tree space can get too large too quickly and slow down the placement for the subsequent samples. Therefore, UShER also provides an option to sort the samples first based on the number of equally parsimonious placements using the `-S` option. 


```
./build/usher -i global_assignments.pb -v <USER_PROVIDED_VCF> -M -S -d output/
```

There are many ways to interpret and visualize the forest of trees produced by multiple placements. One method is to use DensiTree, as shown using an example figure (generated using the [phangorn](https://cran.r-project.org/web/packages/phangorn/) package) below:

<img src="/images/phangorn.png" width="380">

#### Updating multiple input trees

UShER is also fast enough to allow users to update multiple input trees incorporating uncertainty in tree resonstruction, such as multiple bootstrap trees. While we do not provide an explicit option to input multiple trees at once, UShER can be run independently for each input tree and place new samples. We recommend the user to use the [GNU parallel utility](https://www.gnu.org/software/parallel/) to do so in parallel using multiple CPU cores while setting `-T 1` for each UShER task.

## Converting fasta files into a merged VCF for UShER input 

We provide a pipeline for converting SARS-CoV-2 genomic data in fasta format into a merged VCF viable for input to UShER. 

### Alignmenet
If the user possesses SARS-CoV-2 genomic sequences that have been aligned to the reference this step can be skipped. Users can employ multiple alignment using Fast Fourier Transform ([MAFFT](https://mafft.cbrc.jp/alignment/software/)) to perform pairwise alignments between their sequences and the SARS-CoV-2 genomic reference. 

If the user possess separate fasta files for each independent sequence, they must first combine all sequences into one combined fasta file. This can be performed using "cat". For example, for the case of three fasta files named fasta1, fasta2, and fasta3,

```
cat fasta1 fasta2 fasta3 > combined.fa
```

will produce a fasta file named "combined.fa" containing seqences from fasta1, fasta2, and fasta3.

Then, MAFFT can be employed with the following command,

```
mafft --thread 10 --auto --keeplength --addfragments combined.fa reference.fa > myAlignedSequences.fa
```

where reference.fa is a fasta file containing the reference genome and in which case 10 threads are specified. 

### Running faToVcf

Users can use the tool, `faToVcf`, to convert a fasta file containing sequences that have been aligned to a reference into a vcf. For example, if a user possesses a file named "combined.fa" containing sequences that have been aligned to a sequence with fasta header "NC_045512v2", they can run the following command and generate a vcf names "output.vcf":

```
./faToVcf -ref=NC_045512v2 combined.fa output.vcf
```

Users can also mask recommended [problematic sites](https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/14) in SARS-CoV-2 genomic data:

```
./faToVcf -maskSites=problematic_sites_sarsCov2.vcf -ref=NC_045512v2 combined.fa output.vcf
```

The resulting "output.vcf" merged VCF file should be viable for UShER input. Note that in `faToVcf` output, missing data for a particular sample is denoted as "." in the corresponding genotype column.

## matUtils

We are now providing a toolkit, `matUtils`, which can perform a number of tasks related to manipulating and querying the UShER's mutation-annotated tree, such as the generation of the corresponding Newick tree or parsimony-resolved VCF file, masking out mutations, or calculating the number of equally parsimonious placements for a specific set of samples. Full documentation for this toolkit can be found under [src/matUtils](https://github.com/yatisht/usher/blob/master/src/matUtils/README.md).

## Acknowledgement

We thank Jim Kent and the UCSC Genome Browser team for allowing us to download the `faToVcf` utility (from http://hgdownload.soe.ucsc.edu/admin/exe/). Please read the license terms for `faToVcf` here: https://github.com/ucscGenomeBrowser/kent/blob/master/src/LICENSE.

## Reference
**UShER:**
* Yatish Turakhia, Bryan Thornlow, Angie S Hinrichs, Nicola de Maio, Landen Gozashti, Robert Lanfear, David Haussler, and Russ Corbett-Detig, "Ultrafast Sample Placement on Existing Trees (UShER) Empowers Real-Time Phylogenetics for the SARS-CoV-2 Pandemic", bioRxiv [pre-print](https://www.biorxiv.org/content/10.1101/2020.09.26.314971v1) 2020.

**For masking recomendations, please also cite:**
* Landen Gozashti, Conor R. Walker, Nick Goldman, Russell Corbett-Detig and Nicola De Maio, "Issues with SARS-CoV-2 sequencing data: Updated analysis with data from 13th November 2020", Virological 2020 (https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/14
).
