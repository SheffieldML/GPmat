\name{gpsimMapLogLikeGradients}
\Rdversion{1.0}
\alias{gpsimMapLogLikeGradients}
\title{Compute the gradients of the log likelihood of a GPSIMMAP model.}
\description{
  Compute the gradients of the log likelihood of a GPSIMMAP model.
}
\usage{
g <- gpsimMapLogLikeGradients(model)
g <- gpsimMapLogLikeGradients(model, gfdata)
}
\arguments{
  \item{model}{the model for which the log likelihood is computed.}
  \item{gfdata}{the functional gradient of the data portion of the log likelihood}
}
\details{
  \code{g <- gpsimMapLogLikeGradients(model)}
  computes the gradients of the log likelihood of the given Gaussian process for use in a single input motif protein network.

  \code{g <- gpsimMapLogLikeGradients(model, gfdata)}
  computes the gradients of the log likelihood of the given Gaussian process for use in a single input motif protein network.

}
\value{
  \item{g}{the gradients of the parameters of the model.}
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapLogLikelihood}, \link{gpsimMapGradient}, \link{gpsimMapFunctionalLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
