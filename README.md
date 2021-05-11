# Ultrafast Sample Placement on Existing Trees (UShER)

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/yatisht/usher/blob/master/LICENSE

[![License][license-badge]][license-link]
[![Build Status](https://github.com/yatisht/usher/workflows/build/badge.svg)](https://github.com/yatisht/usher/actions)

<img src="/images/usher_logo.png" width="600">

**NEW: We will now be sharing and updating UShER's pre-processed mutation-annotated tree object for public SARS-CoV-2 sequences here: https://hgwdev.gi.ucsc.edu/~angie/UShER_SARS-CoV-2/**

UShER is a program that rapidly places new samples onto an existing phylogeny using maximum parsimony. It is particularly helpful in understanding the relationships of newly sequenced SARS-CoV-2 genomes with each other and with previously sequenced genomes in a global phylogeny. This has emerged as an important challenge during this pandemic for enabling *genomic contact tracing* since the viral phylogeny is already very large (>500K sequences, see https://github.com/roblanf/sarscov2phylo/releases) and is expected to grow by many fold in the coming months. 

UShER is much faster than existing tools with similar functionality and has now also been integrated in the [UCSC SARS-CoV-2 Genome Browser](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace), which does not require UShER installation and usage know-how as described below for SARS-CoV-2 applications. 

Please refer our [wiki](https://usher-wiki.readthedocs.io/) for detailed instructions on installing and using UShER. 

## Acknowledgement

We thank Jim Kent and the UCSC Genome Browser team for allowing us to download the `faToVcf` utility (from http://hgdownload.soe.ucsc.edu/admin/exe/). Please read the license terms for `faToVcf` here: https://github.com/ucscGenomeBrowser/kent/blob/master/src/LICENSE.

## Reference
**UShER:**
* Yatish Turakhia, Bryan Thornlow, Angie S Hinrichs, Nicola de Maio, Landen Gozashti, Robert Lanfear, David Haussler, and Russ Corbett-Detig, "Ultrafast Sample Placement on Existing Trees (UShER) Empowers Real-Time Phylogenetics for the SARS-CoV-2 Pandemic", bioRxiv [pre-print](https://www.biorxiv.org/content/10.1101/2020.09.26.314971v1) 2020.

**For masking recomendations, please also cite:**
* Yatish Turakhia, Nicola De Maio, Bryan Thornlow, Landen Gozashti, Robert Lanfear, Conor R. Walker, Angie S. Hinrichs, Jason D. Fernandes, Rui Borges, Greg Slodkowicz, Lukas Weilguny, David Haussler, Nick Goldman and Russell Corbett-Detig, "Stability of SARS-CoV-2 Phylogenies", PLOS Genetics 2020 (https://doi.org/10.1371/journal.pgen.1009175).
* Landen Gozashti, Conor R. Walker, Robert Lanfear, Nick Goldman, Nicola De Maio and Russell Corbett-Detig, "Issues with SARS-CoV-2 sequencing data: Updated analysis with data from 4 March 2021", Virological 2021 (https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/15).
