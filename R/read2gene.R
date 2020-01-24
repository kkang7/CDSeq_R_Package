#' read2gene outputs the GEP normalized by gene length of the CDSeq-estimated GEPs.
#' \code{read2gene} outputs the gene length normalized CDSeq-estimated GEP.
#' @param read_rate CDSeq-estimated GEP before normalized by gene length.
#' @param gene_effective_length gene effective length which is the gene length minus the read length plus one.
#' @export
#' @return read2gene returns gene length normalized CDSeq-estimated GEPs.

# this function convert read rate to gene rate
# coder: Kai Kang
read2gene <- function(read_rate, gene_effective_length){
  if(nargs()!=2){stop("read2gene function takes 2 inputs")}
  
  g<-nrow(read_rate)
  g1<-length(gene_effective_length)
  if(g!=g1){stop("in function read2gene: two inputs must have the same number of rows")}
  tmp<-read_rate/gene_effective_length
  gene_rate<-t(t(tmp)/colSums(tmp))
  return(gene_rate)
  }