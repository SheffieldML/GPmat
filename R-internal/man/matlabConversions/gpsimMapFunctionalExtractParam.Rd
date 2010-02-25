\name{gpsimMapFunctionalExtractParam}
\Rdversion{1.0}
\alias{gpsimMapFunctionalExtractParam}
\title{Extract the function values from a GPSIMMAP model.}
\description{
  extracts the function values from a structure containing the information about a Gaussian process for single input motif modelling.
}
\usage{
params <- gpsimMapFunctionalExtractParam(model)
}
\arguments{
  \item{model}{the model structure containing the information about the model.}
}
\value{
  \item{params}{a row vector of the function parameters from the model.}
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapFunctionalExpandParam}, \link{modelExtractParam}}.
}
\examples{
## missing
}
\keyword{model}
