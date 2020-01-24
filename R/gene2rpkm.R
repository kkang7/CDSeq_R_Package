#' gene2rpkm outputs the rpkm normalizations of the CDSeq-estimated GEPs.
#' \code{gene2rpkm} outputs the rpkm normalizations of the CDSeq-estimated GEPs.
#' @param gene_rate CDSeq-estimated GEP normalized by gene length.
#' @param gene_effective_length gene effective length which is the gene length minus the read length plus one.
#' @param cell_line_counts RNA-Seq read counts data of cell lines 
#' @return gene2rpkm returns rpkm normalization of the CDSeq-estimated GEPs.

# this function converts the gene rate to rpkm normalization
# coder: Kai Kang

gene2rpkm<-function(gene_rate, gene_effective_length, cell_line_counts){
  if(nargs()!=3){stop("in function gene2rpkm: 3 inputs required")}
  
  #if(!is.matrix(gene_rate)){stop("in function gene2rpkm: input gene_rate has to be a matrix")}
  #if(!is.matrix(cell_line_counts)){stop("in function gene2rpkm: input cell_line_counts has to be a matrix")}
  
  gc<-dim(gene_rate)
  gc2<-dim(cell_line_counts)
  g3<-length(gene_effective_length)
  
  if(gc[1]!=gc2[1] || gc[1]!=g3 || gc2[1]!=g3){stop("in function gene2rpkm: the 3 inputs should have the same number of rows")}
  if(gc[2]!=gc2[2]){stop("in function gene2rpkm: gene_rate and cell_line_counts should have the same number of columns")}
  
  tmp<-cell_line_counts/gene_effective_length
  tmp2<-t(t(tmp)/colSums(cell_line_counts))
  tmp3<-colSums(tmp2)
  rpkm<-t(t(gene_rate)*tmp3)*10^9;
  
  return(rpkm)
}