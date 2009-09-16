\name{gpsimExtractParam}
\Rdversion{1.0}
\alias{gpsimExtractParam}
\title{Extract the parameters of a GPSIM model.}
\description{
  extracts the model parameters from a structure containing the information about a Gaussian process for single input motif modelling.
}
\usage{
params <- gpsimExtractParam(model)
}
\arguments{
  \item{model}{the model structure containing the information about the model.}
}
\value{
  \item{params}{a vector of parameters from the model.}
}
\seealso{
\code{\link{gpsimCreate}, \link{gpsimExpandParam}, \link{modelExtractParam}}.
}
\examples{
## missing
}
\keyword{model}
