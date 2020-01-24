#' CDSeq: A package for complete deconvolution using sequencing data
#' 
#' \code{CDSeq-R-package} takes bulk RNA-seq data as input and simultaneously returns estimates of both cell-type-specific gene expression profiles and sample-specific cell-type proportions.
#' @author Kai Kang, David Huang, \email{kangkai0714@@gmail.com}
#' @references \url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007510}
#' @section Reduce-Recover:
#' CDSeq uses reduce-recovery strategy and CPU parallel computing to speed up the deconvolution.
#' @section Hyperparameter estimation:
#' Estimate hyperparameter for cell-type-specific GEPs (i.e. beta) using reference profile when cell_type_number is scalar.
#' @section Estimating number of cell type:
#' Estimate number of cell types when cell_type_number is a vector of integers.
#' @section Partition on input bulk RNA-seq data:
#' When block_number (number of partition on the bulk RNASeq data) is 1, whole bulk_data will be used. GEP is not from reduce-recovery.
#' @docType package
#' @import Rcpp R.matlab RcppArmadillo MASS foreach doParallel dirmult RcppThread iterators parallel
#' @importFrom Rcpp evalCpp
#' @useDynLib CDSeq
#' @name CDSeq-R-package
NULL 

