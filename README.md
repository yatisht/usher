# Ultrafast Sample Placement on Existing Trees (UShER)

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/yatisht/usher/blob/master/LICENSE

[![License][license-badge]][license-link]
[![Build Status](https://github.com/yatisht/usher/workflows/build/badge.svg)](https://github.com/yatisht/usher/actions)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/usher/README.html)[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=usher)
[![Published in Nature Genetics](https://img.shields.io/badge/Published%20in-Nature%20Genetics-blue.svg)](https://www.nature.com/articles/s41588-021-00862-7)
[![Published in MBE](https://img.shields.io/badge/Published%20in-MBE-blue.svg)](https://doi.org/10.1093/molbev/msab264)
[![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://doi.org/10.1093/bioinformatics/btac401)
[![DOI](https://zenodo.org/badge/296144053.svg)](https://zenodo.org/badge/latestdoi/296144053)



<img src="/images/usher_logo.png" width="600">

**NEW: We will now be sharing and updating UShER's pre-processed mutation-annotated tree object for public SARS-CoV-2 sequences here: http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/. We recommend using https://cov2tree.org/ (developed by [Theo Sanderson](https://github.com/theosanderson/taxodium)) to visualize this tree.**

UShER is now a package consisting of a family of programs for rapid phylogenetic analyses, particularly suitable for the SARS-CoV-2 genomes. 

* **UShER** is a program that rapidly places new samples onto an existing phylogeny using maximum parsimony. It is particularly helpful in understanding the relationships of newly sequenced SARS-CoV-2 genomes with each other and with previously sequenced genomes in a global phylogeny. This has emerged as an important challenge during the COVID-19 pandemic for enabling *genomic contact tracing* since the viral phylogeny is already very large (>2M sequences) and is expected to grow by many fold in the coming months. UShER is much faster than existing tools with similar functionality and has now also been integrated in the [UCSC SARS-CoV-2 Genome Browser](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace), which does not require UShER installation and usage know-how as described below for SARS-CoV-2 applications. If you have sensitive data that cannot be shared over the Internet, consider using [ShUShER](https://shusher.gi.ucsc.edu/), developed by Alex Kramer (https://github.com/amkram/shusher), as an alternative to the Genome Browser. UShER uses the mutation-annotated tree (MAT) data format, which is a phylogenetic tree in which the branches are annotated with the mutations that have been inferred to have occurred on them. 
* **matUtils** is a toolkit for querying, interpreting and manipulating the mutation-annotated trees (MATs). Using matUtils, common operations in SARS-CoV-2 genomic surveillance and contact tracing efforts, including annotating a MAT with new clades, extracting subtrees of the most closely-related samples, or converting the MAT to standard Newick or VCF format can be performed in a matter of seconds to minutes even on a laptop. 
* **matOptimize** is a program to rapidly and effectively optimize a mutation-annotated tree (MAT) for parsimony using subtree pruning and regrafting (SPR) moves within a user-defined radius.
* **RIPPLES** is a program that uses a phylogenomic technique to rapidly and sensitively detect recombinant nodes and their ancestors in a mutation-annotated tree (MAT).  


Please refer to our [wiki](https://usher-wiki.readthedocs.io/) for detailed instructions on installing and using the UShER package. 

## Acknowledgement

We thank Jim Kent and the UCSC Genome Browser team for allowing us to download the `faToVcf` utility (from http://hgdownload.soe.ucsc.edu/admin/exe/). Please read the license terms for `faToVcf` here: https://github.com/ucscGenomeBrowser/kent/blob/master/src/LICENSE.

## References
**UShER:**
* Yatish Turakhia, Bryan Thornlow, Angie S Hinrichs, Nicola de Maio, Landen Gozashti, Robert Lanfear, David Haussler, and Russ Corbett-Detig, "Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic", **Nature Genetics** (2021), [paper](https://t.co/ulGUSRmuWv?amp=1).

**matUtils:**
* Jakob McBroome*, Bryan Thornlow*, Angie S. Hinrichs, Alexander Kramer, Nicola De Maio, Nick Goldman, David Haussler, Russell Corbett-Detig, Yatish Turakhia, "A daily-updated database and tools for comprehensive SARS-CoV-2 mutation-annotated trees", **Molecular Biology and Evolution** (2021), [paper](https://doi.org/10.1093/molbev/msab264).

**RIPPLES:**
* Yatish Turakhia*, Bryan Thornlow*, Angie S. Hinrichs, Jakob McBroome, Nicolas Ayala, Cheng Ye, Nicola De Maio, David Haussler, Russell Corbett-Detig, "Pandemic-Scale Phylogenomics Reveals Elevated Recombination Rates in the SARS-CoV-2 Spike Region", bioRxiv (2021), [preprint](https://www.biorxiv.org/content/10.1101/2021.08.04.455157v1).

**matOptimize**
* Cheng Ye, Bryan Thornlow, Angie Hinrichs, Alexander Kramer, Cade Mirchandani, Devika Torvi, Robert Lanfear, Russell Corbett-Detig, Yatish Turakhia, "matOptimize: A parallel tree optimization method enables online phylogenetics for SARS-CoV-2",  **Bioinformatics** (2021), [paper](https://doi.org/10.1093/bioinformatics/btac401).

**For masking recomendations, please also cite:**
* Yatish Turakhia*, Nicola De Maio*, Bryan Thornlow*, Landen Gozashti, Robert Lanfear, Conor R. Walker, Angie S. Hinrichs, Jason D. Fernandes, Rui Borges, Greg Slodkowicz, Lukas Weilguny, David Haussler, Nick Goldman and Russell Corbett-Detig, "Stability of SARS-CoV-2 Phylogenies", **PLOS Genetics** (2020), [paper](https://doi.org/10.1371/journal.pgen.1009175).
* Landen Gozashti, Conor R. Walker, Robert Lanfear, Nick Goldman, Nicola De Maio and Russell Corbett-Detig, "Issues with SARS-CoV-2 sequencing data: Updated analysis with data from 4 March 2021", Virological 2021 (https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/15).
