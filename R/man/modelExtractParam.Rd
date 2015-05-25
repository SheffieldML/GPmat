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
\alias{gammaPriorExtractParam}
\alias{invgammaPriorExtractParam}
\alias{priorExtractParam}

\title{Extract the parameters of a model.}
\description{
  Extract parameters from the model into a vector of parameters for optimisation.
}
\usage{
  modelExtractParam(model, only.values=TRUE, untransformed.values=FALSE)
}
\arguments{
  \item{model}{the model structure containing the parameters to be extracted.}
  \item{only.values}{include parameter names in the returned vector.}
  \item{untransformed.values}{return actual values, not transformed
    values used by the optimisers.}
}
\value{
  \item{param}{vector of parameters extracted from the model.}
}
\seealso{
  \code{\link{modelExpandParam}}
}
\examples{
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # The probe identifier for TF 'twi'
  twi <- "143396_at"
  # The probe identifier for the target gene
  targetProbe <- "152715_at"

  # Create the model, but do not optimise
  model <- GPLearn(drosophila_gpsim_fragment,
                   TF=twi, targets=targetProbe,
                   useGpdisim=TRUE, quiet=TRUE,
                   dontOptimise=TRUE)

  # Get the initial parameter values
  params <- modelExtractParam(model, only.values=FALSE)
}
\keyword{model}
