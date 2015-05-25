\name{drosophila_gpsim_fragment}
\docType{data}
\alias{drosophila_gpsim_fragment}
\title{Fragment of 12 time point Drosophila embryonic development
  microarray gene expression time series}
\description{
  Four genes from the 12 time point Drosophila embryonic development
  Affymetrix microarray gene expression data set by Tomancak et
  al. (2002).

  The data has been processed using \code{\link{mmgmos}} and
  \code{\link{processData}}.
}
\usage{data(drosophila_gpsim_fragment)}
\format{An \code{\linkS4class{ExpressionTimeSeries}} object with
  3 repeats of the 12 time points for 4 probes.}
\source{\url{ftp://ftp.fruitfly.org/pub/embryo_tc_array_data/}}
\references{
  Tomancak, P et al.
  Systematic determination of patterns of gene expression during
  Drosophila embryogenesis.
  \emph{Genome Biol} 3:RESEARCH0088, 2002.
}
\keyword{datasets}
