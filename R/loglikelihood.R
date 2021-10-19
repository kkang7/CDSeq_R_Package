#' loglikelihood computes the log likelihood of the CDSeq model.
#' \code{loglikelihood} outputs the value of log likelihood
#' @param estProp CDSeq-estimated cell type proportions.
#' @param estGEP CDSeq-estimated cell-type-specific GEPs.
#' @param mydata input bulk RNA-seq data.
#' @param alpha hyperparameter for cell type proportion estimation.
#' @param beta hyperparameter for cell-type-specific GEP estimation.
#' @return loglikelihood returns log likelihood values.

# coder: Kai Kang

loglikelihood <- function(estProp, estGEP, mydata, alpha, beta){
  #if(class(alpha)!="numeric" || length(alpha)>1){stop("CDSeq: in function loglikelihood, alpha is a positive scalar")}
  #if(class(beta)!="numeric" || length(beta)>1){stop("beta is a positive scalar")}
  if(class(alpha)!="numeric"){stop("CDSeq: in function loglikelihood, alpha should be a numeric value")}
  #if(class(beta)!="numeric" && sum(!("matrix" %in% class(beta)))){stop("CDSeq: in function loglikelihood, beta should be a numeric value")}
  mu<-1e8
  delta <- 1e-15 # to avoid Inf value of log
  # value <- sum(log(mu*estGEP)*(t(beta) - 1)) + sum(log(estProp)*(alpha - 1)) + sum((log(mu*estGEP %*% estProp  )) * mydata)
  # if(is.infinite(value)){
  #   value <- sum(log(mu* (estGEP+delta))*(t(beta) - 1)) + sum(log(estProp)*(alpha - 1)) + sum((log(mu*(estGEP+delta) %*% estProp  )) * mydata)
  # }
  
  value <- sum((log(mu*(estGEP+delta) %*% estProp  )) * mydata)
  return(value)
}