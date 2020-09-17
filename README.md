# Ultrafast Sample Placement on Existing Trees (UShER)

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/yatisht/usher/blob/master/LICENSE

[![License][license-badge]][license-link]
[![Build Status](https://github.com/yatisht/usher/workflows/build/badge.svg)](https://github.com/yatisht/usher/actions)

<img src="/images/usher_logo.png" width="600">

UShER is a program that rapidly places new samples onto an existing phylogeny using maximum parsimony. It is particularly helpful in understanding the relationships of newly sequenced SARS-CoV-2 genomes with each other and with previously sequenced genomes in a global phylogeny. This has emerged as an important challenge during this pandemic for enabling *genomic contact tracing* since the viral phylogeny is already very large (>50K sequences, see https://github.com/roblanf/sarscov2phylo/releases) and is expected to grow by many fold in the coming months. 

UShER is much faster than existing tools with similar functionality and has now also been integrated in the [UCSC SARS-CoV-2 Genome Browser](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace).

## Installing UShER

First, clone the UShER repository using the commands below.
```
git clone https://github.com/yatisht/strain_phylogenetics.git
cd strain_phylogenetics
```

Next, you may install UShER either using [Docker](https://www.docker.com/) (**recommended**) or using one of the installation scripts that we provide for **MacOS (10.14 and above)**, **Ubuntu (18.04 and above)** and **CentOS (7 and above)**. See instructions below.

### Using Docker

Perhaps the simplest method to install this package irrespective of the platform (Windows, Linux or MacOS) is using [Docker](https://www.docker.com/) (see installation instructions [here](https://docs.docker.com/get-docker/)). Once Docker is installed (which takes a few minutes typically), use the following command to build and install UShER.

```
docker build -t usher .
```
Once installed, you may start a new bash in a Ubuntu 18.04 container using the Docker command below, following which the remaining commands on this page should work smoothly.
```
docker run -t -i usher /bin/bash
```

### Using installation scripts

If you don't have Docker or prefer not to use it for some reason, you may use one of our installation scripts below to install UShER.

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

## Using UShER

### Pre-processing global phylogeny

### Placing new samples

### Uncertainty in placing new samples

* UShER is a program that rapidly places new samples onto an existing phylogeny using maximum parsimony. It is particularly helpful in understanding the relationships of newly sequenced SARS-CoV-2 genomes with each other and with previously sequenced genomes in an existing phylogeny. UShER prep-processes the existing phylogeny (pruned-sumtree-for-cog.nh in the example below) and its variation (pruned-sumtree-for-cog.vcf in the example below), computes the parsimonious assignments of each variation and stores the results in a compact [protobuf](https://developers.google.com/protocol-buffers) file (pruned-sumtree-for-cog.assignments.pb in the example below). 
```
    $ ./build/usher --tree tree/pruned-sumtree-for-cog.nh --vcf vcf/pruned-sumtree-for-cog.vcf --threads 4 --save-assignments pruned-sumtree-for-cog.assignments.pb 
```
* Once the pre-processing is complete, new sequences whose variants are called in a VCF file (missing.vcf in the example below) can be added to the existing phylogeny using the command below:
```
    $ ./build/usher --load-assignments pruned-sumtree-for-cog.assignments.pb --vcf vcf/missing.vcf  --threads 4
```
* UShER is much faster than existing tools with similar functionality and has now been integrated in the [UCSC SARS-CoV-2 Genome Browser](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace).

## Reference
* Coming soon
