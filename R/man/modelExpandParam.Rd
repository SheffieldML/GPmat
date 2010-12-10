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
\alias{modelUpdateProcesses}
\alias{cgpdisimUpdateProcesses}
\alias{gpdisimUpdateProcesses}
\alias{cgpsimUpdateProcesses}
\alias{gpsimUpdateProcesses}

\title{Update a model structure with new parameters or update the
  posterior processes.}
\description{
  Update a model structure or component with new parameters,
  or update the posterior processes.
}
\usage{
model <- modelExpandParam(model, param, untransformed.values=FALSE)
model <- modelUpdateProcesses(model)
}
\arguments{
  \item{model}{the model structure to be updated.}
  \item{param}{vector of parameters.}
  \item{untransformed.values}{vector contains actual values, not
    transformed values used by the optimisers.}
}
\value{
  \item{model}{updated model structure.}
}
\details{
  \code{model <- modelExpandParam(model, param)}
  returns a model structure filled with the parameters in the given
  vector. This is used as a helper function to enable parameters to be
  optimised in, for example, the optimisation functions.
  
  \code{model <- modelUpdateProcesses(model)}
  updates posterior processes of the given model.
}
\seealso{
\code{
  \link{modelExtractParam}}
}
\keyword{model}
