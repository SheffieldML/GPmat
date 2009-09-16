\name{gpsimExpandParam}
\Rdversion{1.0}
\alias{gpsimExpandParam}
\title{Expand the given parameters into a GPSIM structure.}
\description{
  takes the given vector of parameters and places them in the model structure, it then updates any stored representations that are dependent on those parameters, for example kernel matrices etc..
}
\usage{
model <- gpsimExpandParam(model, params)
}
\arguments{
  \item{model}{the model structure for which parameters are to be updated.}
  \item{params}{a vector of parameters for placing in the model structure.}
}
\value{
  \item{model}{a returned model structure containing the updated parameters. }
}
\seealso{
\code{\link{gpsimCreate}, \link{gpsimExtractParam}, \link{modelExtractParam}, \link{gpsimUpdateKernels}}.
}
\examples{
## missing
}
\keyword{model}
