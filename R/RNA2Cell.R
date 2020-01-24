#' RNA proportion to cell proportion
#' \code{RNA2Cell} converts RNA proportion to cell proportion
#' @param eta numeric vector represents the different amounts of RNA produced by different cell types
#' @param rnaprop sample-specific cell-type RNA proportion
#' @return RNA2Cell returns sample-specific cell-type proportion
# this function converts RNA proportion to cell proportion
# coder: Kai Kang

RNA2Cell<-function(eta,rnaprop){
  if(nargs()!=2){stop("in RNA2Cell: function takes 2 inputs")}
  
  if(!is.null(dim(eta))){stop("in RNA2Cell: eta should be a vector, NOT a matrix")}
  nc<-length(eta)
  nc2<-nrow(rnaprop)
  if(is.null(nc2)){stop("in RNA2Cell: rnaprop should be a matrix")}
  if(nc!=nc2){stop("in RNA2Cell: length(eta) should be equal to nrow(rnaprop)")}
  
  tmp<-rnaprop/eta
  cellprop<-t(t(tmp)/colSums(tmp))
  return(cellprop)
}