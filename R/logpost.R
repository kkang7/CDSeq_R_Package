#' logpost computes the log posterior of the CDSeq model.
#' \code{logpost} outputs the value of log posterior.
#' @param estProp CDSeq-estimated cell type proportions.
#' @param estGEP CDSeq-estimated cell-type-specific GEPs.
#' @param mydata input bulk RNA-seq data.
#' @param alpha hyperparamter for cell type proportion estimation.
#' @param beta hyperparameter for cell-type-specific GEP estimation.
#' @return logpost returns log posterior values.

## simplified log posterior function with estProp and estGEP matrices as input
# coder: Kai Kang

logpost <- function(estProp, estGEP, mydata, alpha, beta){
  #if(class(alpha)!="numeric" || length(alpha)>1){stop("CDSeq: in function logpost, alpha is a positive scalar")}
  #if(class(beta)!="numeric" || length(beta)>1){stop("beta is a positive scalar")}
  if(class(alpha)!="numeric"){stop("CDSeq: in function logpost, alpha should be a numeric value")}
  if(class(beta)!="numeric"){stop("CDSeq: in function logpost, beta should be a numeric value")}
  mu<-1e8
  value <- sum(log(mu*estGEP)*(beta - 1)) + sum(log(estProp)*(alpha - 1)) + sum((log(mu*estGEP %*% estProp  )) * mydata)
  return(value)
}