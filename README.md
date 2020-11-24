# QuRIE-seq_manuscript

## Abstract

Current high-throughput single-cell multi-omics methods cannot concurrently map changes in (phospho)protein levels and the associated gene expression profiles. We present QuRIE-seq (Quantification of RNA and Intracellular Epitopes by sequencing) and use multi-factor omics analysis (MOFA+) to map signal transduction over multiple timescales. We demonstrate that QuRIE-seq can trace the activation of the B-cell receptor pathway at the minute and hour time-scale and provide insight into the mechanism of action of an inhibitory drug Ibrutinib.

## In this repository

### Folder structure

Build according to a [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr

To run the code: clone the repository. Download the raw count tables from GEO #### and save in data/raw/ folder. 

First run the QC.Rmd file. After that either MOFAaIg or MOFAibru Rmd files can be used. 

### Github page 

The docs folder contains the actual github page. 

Three pages contain code to process, analyse and create figures:

* [QC](QC.html) processes count tables (available at GEO####): filter high quality cells; normalize and scale counts.   
* [MOFA aIg](MOFAaIg.html) computes a model of cells stimulates with aIg molecule.   
* [MOFA ibru](MOFAibru.html) computes a model of cells stimulates with aIg molecule with or without prescence of ibrutinib.  


