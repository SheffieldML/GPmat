\name{gpsimMapFunctionalExpandParam}
\Rdversion{1.0}
\alias{gpsimMapFunctionalExpandParam}
\title{Expand the given function values into a GPSIMMAP structure.}
\description{
  takes the given vector of function values and places them in the model structure. It then updates any stored representations that are dependent on those function values, for example the Hessian, etc..
}
\usage{
model <- gpsimMapFunctionalExpandParam(model, vals)
}
\arguments{
  \item{model}{the model structure for which parameters are to be updated.}
  \item{vals}{a vector of function values for placing in the model structure.}
}
\value{
  \item{model}{a returned model structure containing the updated function values. }
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapFunctionalExtractParam}, \link{modelExtractParam}, \link{gpsimMapUpdateKernels}}.
}
\examples{
## missing
}
\keyword{model}
