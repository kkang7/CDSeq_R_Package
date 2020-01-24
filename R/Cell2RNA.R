#' Cell proportion to RNA proportion
#' \code{Cell2RNA} converts Cell proportion to RNA proportion
#' @param eta numeric vector represents the different amounts of RNA produced by different cell types
#' @param cellprop sample-specific cell-type proportion
#' @return Cell2RNA returns sample-specific cell-type RNA proportion

# this function converts RNA proportion to cell proportion
# coder: Kai Kang

Cell2RNA<-function(eta,cellprop){
  if(nargs()!=2){stop("in Cell2RNA: function takes 2 inputs")}
  
  if(!is.null(dim(eta))){stop("in Cell2RNA: eta should be a vector, NOT a matrix")}
  nc<-length(eta)
  nc2<-nrow(cellprop)
  if(is.null(nc2)){stop("in Cell2RNA: cellprop should be a matrix")}
  if(nc!=nc2){stop("in Cell2RNA: length(eta) should be equal to nrow(cellprop)")}
  
  tmp<-cellprop*eta
  rnaprop<-t(t(tmp)/colSums(tmp))
  return(rnaprop)
}