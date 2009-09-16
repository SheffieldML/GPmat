\name{gpsimMapFunctionalLikeGrad2}
\Rdversion{1.0}
\alias{gpsimMapFunctionalLikeGrad2}
\title{Compute the functional gradient for GPSIMMAP.}
\description{
  computes the functional gradient of the log likelihood for use in the MAP approximation to the GPSIM posterior solution.
}
\usage{
p <- gpsimMapFunctionalLikeGrad2(model)
}
\arguments{
  \item{model}{the model for which the gradient is to be computed.}
}
\value{
  \item{p}{the gradient of the log likelihood with respect to the points of the function. }
}
\seealso{
\code{\link{gpsimMapLikeGradientImplicit}}.
}
\examples{
## missing
}
\keyword{model}
