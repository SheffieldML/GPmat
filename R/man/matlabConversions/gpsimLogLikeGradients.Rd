\name{gpsimLogLikeGradients}
\Rdversion{1.0}
\alias{gpsimLogLikeGradients}
\title{Compute the gradients of the log likelihood of a GPSIM model.}
\description{
  computes the gradients of the log likelihood of the given Gaussian process for use in a single input motif protein network.
}
\usage{
g <- gpsimLogLikeGradients(model)
}
\arguments{
  \item{model}{the model for which the log likelihood is computed.}
}
\value{
  \item{g}{the gradients of the parameters of the model. }
}
\seealso{
\code{\link{gpsimCreate}, \link{gpsimLogLikelihood}, \link{gpsimGradient}}.
}
\examples{
## missing
}
\keyword{model}
