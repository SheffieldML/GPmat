\name{gpsimMapLogLikelihood}
\Rdversion{1.0}
\alias{gpsimMapLogLikelihood}
\title{Compute the log likelihood of a GPSIMMAP model.}
\description{
  computes the log likelihood of the given Gaussian process for use in a single input motif protein network.
}
\usage{
ll <- gpsimMapLogLikelihood(model)
}
\arguments{
  \item{model}{the model for which the log likelihood is computed.}
}
\value{
  \item{ll}{the log likelihood of the data set. }
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapLogLikeGradient}, \link{gpsimMapObjective}}.
}
\examples{
## missing
}
\keyword{model}
