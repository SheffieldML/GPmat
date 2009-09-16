\name{gpsimMapExpandParam}
\Rdversion{1.0}
\alias{gpsimMapExpandParam}
\title{Expand the given parameters into a GPSIMMAP structure.}
\description{
  takes the given vector of parameters and places them in the model structure, it then updates any stored representations that are dependent on those parameters, for example kernel matrices etc..
}
\usage{
model <- gpsimMapExpandParam(model, params)
}
\arguments{
  \item{model}{the model structure for which parameters are to be updated.}
  \item{params}{a vector of parameters for placing in the model structure.}
}
\value{
  \item{model}{a returned model structure containing the updated parameters. }
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapExtractParam}, \link{modelExtractParam}, \link{gpsimMapUpdateKernels}}.
}
\examples{
## missing
}
\keyword{model}
