# CDSeq R Package
CDSeq is a complete deconvolution method for dissecting bulk RNA-Seq data. The input of CDSeq is, ideally, bulk RNA-Seq read counts (for example HTSeq read counts), and CDSeq will estimate, simultaneously, the cell-type-specific gene expression profiles and the sample-specific cell-type proportions, no reference of pure cell line GEPs or scRNAseq reference is needed for running CDSeq. 

For example, if you have a bulk RNA-Seq data, $A_{G\times M}$, which is a G by M matrix. G denotes the number of genes and M is the sample size, then CDSeq will output $B_{G\times T}$ and $C_{T\times M}$, where T is the number of cell types, $B_{G\times T}$ is the estiamte of cell-type-specific GEPs and $C_{T\times M}$ is the estiamte of sample-specific cell-type proportions. 

Importantly, you can ask CDSeq to estimate the number of cell types, i.e. T, by providing a vector of possible integer values for T. For example, if the user input for T is a vector, i.e. $T=\{2,3,4,5,6\}$, then CDSeq will estimate the most likely number for T.    

 
## Versioning 
I follow the versioning scheme: *major*.*minor*.*patch*. In CDSeq package, *major*, *minor* and *patch* basically indicate the following types of changes:
1. *major*: changes on the underlying probabilistic modeling 
2. *minor*: changes on the numerical procedures or coding structure
3. *patch*: others


## Version 1.0.8
Bug fix. Add vignette.

## Version 1.0.7
Bug fix. Add more features to cellTypeAssignSCRNA function.

## Version 1.0.6
Add cellTypeAssignSCRNA function. This function can be used for cell type assignment using user provided single cell RNAseq data (with annotations).


## Version 1.0.5
Add cellTypeAssignMarkerGenes function. This function can be used for cell type assignment using user provided marker gene list.



## Installation

install_github("kkang7/CDSeq_R_Package")

## scRNAseq data for CDSeq-estimated cell type annotations
We collected some public available scRNAseq data from https://www.nature.com/articles/s41586-020-2157-4 and put them in hdf5 format. The data are available upon request. 

## Known issue about MacOS installation
It is possible for Mac users to run into some errors when install from source due to problems of Rcpp compiler tools. Follow the instruction here may help: https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/
 

## Contact 
email: kangkai0714@gmail.com
