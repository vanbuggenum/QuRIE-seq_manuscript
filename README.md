# QuRIE-seq_manuscript

## Abstract

Current high-throughput single-cell multi-omics methods cannot concurrently map changes in (phospho)protein levels and the associated gene expression profiles. We present QuRIE-seq (Quantification of RNA and Intracellular Epitopes by sequencing) and use multi-factor omics analysis (MOFA+) to map signal transduction over multiple timescales. We demonstrate that QuRIE-seq can trace the activation of the B-cell receptor pathway at the minute and hour time-scale and provide insight into the mechanism of action of an inhibitory drug Ibrutinib.

## In this repository

### Folder structure

Build using [workflowr](https://github.com/jdblischak/workflowr) project.

A nice html-page that is generated from the Rmd files in the analysis folder can be accessed [here](https://vanbuggenum.github.io/QuRIE-seq_manuscript/index.html)

To run the Rmd files: 1) clone the repository. 2) Download the raw count tables from [GSE162461](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162461) and save in data/raw/ folder. 3) First run the QC.Rmd file. After that either MOFAaIg or MOFAibru Rmd files can be used. 

### Github page 

The docs folder contains the actual github page. 

Three pages contain code to process, analyse and create figures:

* [QC](https://vanbuggenum.github.io/QuRIE-seq_manuscript/QC.html) processes count tables (available at [GSE162461](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162461): filter high quality cells; normalize and scale counts.   
* [MOFA aIg](https://vanbuggenum.github.io/QuRIE-seq_manuscript/MOFAaIg.html) computes a model of cells stimulates with aIg molecule.   
* [MOFA ibru](https://vanbuggenum.github.io/QuRIE-seq_manuscript/MOFAibru.html) computes a model of cells stimulates with aIg molecule with or without prescence of ibrutinib.  


