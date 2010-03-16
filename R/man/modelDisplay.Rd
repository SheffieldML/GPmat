\name{modelDisplay}
\Rdversion{1.0}
\alias{modelDisplay}
\alias{gpdisimDisplay}
\alias{gpsimDisplay}
\alias{cmpndKernDisplay}
\alias{disimKernDisplay}
\alias{kernDisplay}
%\alias{mlpKernDisplay}
\alias{multiKernDisplay}
\alias{rbfKernDisplay}
\alias{simKernDisplay}
%\alias{translateKernDisplay}
\alias{whiteKernDisplay}
\title{Display a model.}
\description{
  displays the parameters of the model/kernel and the model/kernel type to
  the console.
}
\usage{
model <- modelDisplay(model)
}
\arguments{
  \item{model}{the model/kernel structure to be displayed.}
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
  #modelDisplay(model)
}
\keyword{model}
