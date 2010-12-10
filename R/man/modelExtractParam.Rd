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
  param <- modelExtractParam(model, only.values=TRUE, untransformed.values=FALSE)
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
\keyword{model}
