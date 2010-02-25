\name{gpsimMapGradients}
\Rdversion{1.0}
\alias{gpsimMapGradients}
\title{Compute the gradients of the log likelihood of a GPSIMMAP model.}
\description{
  computes the gradients of the log likelihood of the given Gaussian process for use in a single input motif protein network.
}
\usage{
g <- gpsimMapGradients(model)
}
\arguments{
  \item{model}{the model for which the log likelihood is computed.}
}
\value{
  \item{g}{the gradients of the parameters of the model. }
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapLogLikelihood}, \link{gpsimMapGradient}, \link{gpsimMapFunctionalLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
