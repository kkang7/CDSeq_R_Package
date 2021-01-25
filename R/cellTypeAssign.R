#' Assign cell types using correlation matrix computed using cell-type-specific GEPs and reference GEPs.
#' \code{cellTypeAssign} assigns CDSeq-identified cell types to reference profile.
#' @param corMat correlation matrix between CDSeq-estimated GEPs and reference GEPs.
#' @param threshold only the correlations that are above threshold will be considered.
#' @return cellTypeAssign returns a vector of cell type assignment to the reference profile. 


# cell type assignment using correlation matrix
# coder: Kai Kang
# last updated: 5/3/2019

cellTypeAssign<-function(corMat,threshold=0.8){
  
  # corMat is the correlation matrix. 
  # dimension of corMat is estimated cell type GEPs by reference GEP profile
  n<-nrow(corMat)
  
  cellTypeAssign_result<-rep(0,n)
  for (i in 1:n) {
    if(max(corMat[i,])>threshold){cellTypeAssign_result[i]<-which.max(corMat[i,])}
  }
  return(cellTypeAssign_result)
}