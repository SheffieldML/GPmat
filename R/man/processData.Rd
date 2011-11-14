\name{processData}
\Rdversion{1.0}
\alias{processData}
\alias{processRawData}
\title{Processing expression time series}
\description{
  \code{processData} further processes time series data preprocessed
  by \code{puma} or \code{lumi}.
  
  \code{processRawData} performs similar processing for other data.

  Both functions return \code{\linkS4class{ExpressionTimeSeries}}
  objects that can be used as input for the functions
  \code{\link{GPLearn}} and \code{\link{GPRankTargets}}.
}
\usage{
  processData(data, times = NULL, experiments = NULL,
    do.normalisation = TRUE)
  processRawData(rawData, times, experiments = NULL,
    is.logged = TRUE, do.normalisation = ifelse(is.logged, TRUE, FALSE))
}
\arguments{
  \item{data}{The preprocessed data from \code{mmgMOS} in the
    \code{puma} package (an \code{exprReslt} object) or
    from the \code{lumi} package (a \code{LumiBatch} object).}
  \item{rawData}{Raw data matrix to be used. Each row corresponds to
    a gene and each column to a data point.}
  \item{times}{Observation times of each data point.
    If unspecified or NULL, \code{processData} attempts to infer this
    from phenoData(data) field containing 'time' in the name.}
  \item{experiments}{The replicate structure of the data indicating
    which expression data points arise from which experiments.
    This should be an array in integers from 1 to N with length
    equal to the number of data points.
    By default all the data points are assumed to be from same
    replicate.}
  \item{is.logged}{Indicates whether the expression values are on
    log scale or not. Normalisation of non-logged data is unsupported.}
  \item{do.normalisation}{Indicates whether to perform the
    normalisation.}
}
\details{
  The expression data (and percentiles, if available) are normalized
  by equalising the mean of log-expression in each time points.
  In \code{processData}, a normal
  distribution is then fitted into the data with distfit.
}
\value{
  An \code{\linkS4class{ExpressionTimeSeries}}
  object containing all provided information.
}
\author{Antti Honkela, Jonatan Ropponen}
\seealso{
  \code{\link{GPLearn}, \link{GPRankTargets}}.
}
\examples{
  ## Load a mmgmos preprocessed fragment of the Drosophila developmental
  ## time series
  data(drosophila_mmgmos_fragment)

  ## Process the data (3 experiments containing 12 time points each)
  drosophila_gpsim_fragment <- processData(drosophila_mmgmos_fragment,
    experiments=rep(1:3, each=12))
}
\keyword{model}
