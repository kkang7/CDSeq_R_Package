% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logpost.R
\name{logpost}
\alias{logpost}
\title{logpost computes the log posterior of the CDSeq model.
\code{logpost} outputs the value of log posterior.}
\usage{
logpost(estProp, estGEP, mydata, alpha, beta)
}
\arguments{
\item{estProp}{CDSeq-estimated cell type proportions.}

\item{estGEP}{CDSeq-estimated cell-type-specific GEPs.}

\item{mydata}{input bulk RNA-seq data.}

\item{alpha}{hyperparameter for cell type proportion estimation.}

\item{beta}{hyperparameter for cell-type-specific GEP estimation.}
}
\value{
logpost returns log posterior values.
}
\description{
logpost computes the log posterior of the CDSeq model.
\code{logpost} outputs the value of log posterior.
}
