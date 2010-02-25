\name{gpsimMapFunctionalLogLikeGradients}
\Rdversion{1.0}
\alias{gpsimMapFunctionalLogLikeGradients}
\title{Compute the functional gradient for GPSIMMAP.}
\description{
  computes the functional gradient of the log likelihood for use in the MAP approximation to the GPSIM posterior solution.
}
\usage{
g <- gpsimMapFunctionalLogLikeGradients(model)
}
\arguments{
  \item{model}{the model for which the gradient is to be computed.}
}
\value{
  \item{g}{the gradient of the log likelihood with respect to the points of the function. }
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapUpdateYpred}}.
}
\examples{
## missing
}
\keyword{model}
