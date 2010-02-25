\name{gpsimMapExtractParam}
\Rdversion{1.0}
\alias{gpsimMapExtractParam}
\title{Extract the parameters of a GPSIMMAP model.}
\description{
  extracts the model parameters from a structure containing the information about a Gaussian process for single input motif modelling.
}
\usage{
params <- gpsimMapExtractParam(model)
}
\arguments{
  \item{model}{the model structure containing the information about the model.}
}
\value{
  \item{params}{a vector of parameters from the model.}
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapExpandParam}, \link{modelExtractParam}}.
}
\examples{
## missing
}
\keyword{model}
