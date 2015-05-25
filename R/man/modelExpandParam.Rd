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
\alias{gammaPriorExpandParam}
\alias{invgammaPriorExpandParam}
\alias{priorExpandParam}

\title{Update a model structure with new parameters or update the
  posterior processes.}
\description{
  Update a model structure or component with new parameters,
  or update the posterior processes.
}
\usage{
  modelExpandParam(model, params)
  modelUpdateProcesses(model, predt=NULL)
}
\arguments{
  \item{model}{the model structure to be updated.}
  \item{params}{vector of parameters.}
  \item{predt}{(optional) a vector of times to infer the posterior at.
    By default this is 100 points spanning the time range of the
    data.}
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
\code{\link{GPLearn}, \link{modelExtractParam}}
}
\examples{\dontrun{
  # Learn the model
  model <- GPLearn(...)
  params <- modelExtractParam(model, only.values=TRUE)
  params[1] <- 0
  new_model <- modelExpandParam(model, params)
  new_model <- modelUpdateProcesses(new_model)
}}
\keyword{model}
