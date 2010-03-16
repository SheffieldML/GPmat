\name{modelExpandParam}
\Rdversion{1.0}
\alias{modelExpandParam}
\alias{cgpdisimExpandParam}
\alias{cgpsimExpandParam}
\alias{gpdisimExpandParam}
\alias{gpsimExpandParam}
\alias{cmpndKernExpandParam}
\alias{disimKernExpandParam}
\alias{kernExpandParam}
\alias{mlpKernExpandParam}
\alias{multiKernExpandParam}
\alias{rbfKernExpandParam}
\alias{simKernExpandParam}
\alias{translateKernExpandParam}
\alias{whiteKernExpandParam}
\title{Update a model structure with parameters.}
\description{
  returns a model structure filled with the parameters in the given vector. This is used as a helper function to enable parameters to be optimised in, for example, the NETLAB optimisation functions.
}
\usage{
model <- modelExpandParam(model, param)
}
\arguments{
  \item{model}{the model structure in which the parameters are to be placed.}
  \item{param}{vector of parameters which are to be placed in the model structure.}
}
\value{
  \item{model}{model structure with the given parameters in the relevant locations.}
}
\seealso{
\code{\link{modelExtractParam}}
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
  #params <- modelExtractParam(model, only.values=TRUE)
  #params[1] <- 0
  #new_model <- modelExpandParam(model, params)
}
\keyword{model}
