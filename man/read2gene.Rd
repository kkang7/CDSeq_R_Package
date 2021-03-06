% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read2gene.R
\name{read2gene}
\alias{read2gene}
\title{read2gene outputs the GEP normalized by gene length of the CDSeq-estimated GEPs.
\code{read2gene} outputs the gene length normalized CDSeq-estimated GEP.}
\usage{
read2gene(read_rate, gene_effective_length)
}
\arguments{
\item{read_rate}{CDSeq-estimated GEP before normalized by gene length.}

\item{gene_effective_length}{gene effective length which is the gene length minus the read length plus one.}
}
\value{
read2gene returns gene length normalized CDSeq-estimated GEPs.
}
\description{
read2gene outputs the GEP normalized by gene length of the CDSeq-estimated GEPs.
\code{read2gene} outputs the gene length normalized CDSeq-estimated GEP.
}
