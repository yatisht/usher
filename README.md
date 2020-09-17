# Ultrafast Sample Placement on Existing Trees (UShER)

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/yatisht/usher/blob/master/LICENSE

[![License][license-badge]][license-link]
[![Build Status](https://github.com/yatisht/usher/workflows/build/badge.svg)](https://github.com/yatisht/usher/actions)

<img src="/images/usher_logo.png" width="600">

### Install prerequisites and build programs
* Perhaps the simplest build process for this package on any platform (Windows, Linux or MacOS) is using [Docker](https://www.docker.com/) (see installation instructions [here](https://docs.docker.com/get-docker/)). Once Docker is installed (which takes a few minutes typically), the following two commands would build the programs and start a new bash in an Ubuntu 18.04 container following which the remaning commands on this page should work smoothly.
```
    $ docker build -t usher .
    $ docker run -t -i usher /bin/bash
```
* UShER is a program that rapidly places new samples onto an existing phylogeny using maximum parsimony. It is particularly helpful in understanding the relationships of newly sequenced SARS-CoV-2 genomes with each other and with previously sequences genomes in an existing phylogeny. UShER prep-processes the existing phylogeny (pruned-sumtree-for-cog.nh in the example below) and its variation (pruned-sumtree-for-cog.vcf in the example below), computes the parsimonious assignments of each variation and stores the results in a compact [protobuf](https://developers.google.com/protocol-buffers) file (pruned-sumtree-for-cog.assignments.pb in the example below). 
```
    $ ./build/usher --tree tree/pruned-sumtree-for-cog.nh --vcf vcf/pruned-sumtree-for-cog.vcf --threads 4 --save-assignments pruned-sumtree-for-cog.assignments.pb 
```
* Once the pre-processing is complete, new sequences whose variants are called in a VCF file (missing.vcf in the example below) can be added to the existing phylogeny using the command below:
```
    $ ./build/usher --load-assignments pruned-sumtree-for-cog.assignments.pb --vcf vcf/missing.vcf  --threads 4
```
* UShER is much faster than existing tools with similar functionality and has now been integrated in the [UCSC SARS-CoV-2 Genome Browser](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace).

