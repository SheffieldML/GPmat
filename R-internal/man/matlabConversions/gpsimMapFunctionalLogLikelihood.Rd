\name{gpsimMapFunctionalLogLikelihood}
\Rdversion{1.0}
\alias{gpsimMapFunctionalLogLikelihood}
\title{Compute the log likelihood of a GPSIMMAP model.}
\description{
  computes the log likelihood of the given Gaussian process for use in a single input motif protein network. Acts as a wrapper for gpsimMapLogLikelihood.
}
\usage{
ll <- gpsimMapFunctionalLogLikelihood(model)
}
\arguments{
  \item{model}{the model for which the log likelihood is computed.}
}
\value{
  \item{ll}{the log likelihood of the data set. }
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapLogLikelihood}, \link{gpsimMapFunctionalLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
