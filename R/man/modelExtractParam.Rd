\name{modelExtractParam}
\Rdversion{1.0}
\alias{modelExtractParam}
\alias{cgpdisimExtractParam}
\alias{cgpsimExtractParam}
\alias{gpdisimExtractParam}
\alias{gpsimExtractParam}
\alias{cmpndKernExtractParam}
\alias{disimKernExtractParam}
\alias{kernExtractParam}
\alias{mlpKernExtractParam}
\alias{multiKernExtractParam}
\alias{rbfKernExtractParam}
\alias{simKernExtractParam}
\alias{translateKernExtractParam}
\alias{whiteKernExtractParam}
\title{Extract the parameters of a model.}
\description{
  Extract parameters from the model into a vector of parameters for optimisation.
}
\usage{
  param <- modelExtractParam(model, only.values=TRUE)
}
\arguments{
  \item{model}{the model structure containing the parameters to be extracted.}
  \item{only.values}{include parameter names in the returned vector.}
}
\value{
  \item{param}{vector of parameters extracted from the model.}
}
\seealso{
  \code{\link{modelExpandParam}}
}
\examples{
  ## Commented out to speed up package check
  ## Load a mmgmos preprocessed fragment of the Drosophila developmental
  ## time series
  #data(drosophila_gpsim_fragment)

  ## Get the target probe names
  #library(annotate)
  #aliasMapping <- getAnnMap("ALIAS2PROBE",
  #                  annotation(drosophila_gpsim_fragment))
  #twi <- get('twi', env=aliasMapping)
  #fbgnMapping <- getAnnMap("FLYBASE2PROBE",
  #                 annotation(drosophila_gpsim_fragment))
  #targetProbe <- get('FBgn0035257', env=fbgnMapping)

  ## Learn the model
  #model <- GPLearn(drosophila_gpsim_fragment,
  #                 TF=twi, targets=targetProbe,
  #                 useGpdisim=TRUE, quiet=TRUE)
  #params <- modelExtractParam(model, only.values=FALSE)
}
\keyword{model}
