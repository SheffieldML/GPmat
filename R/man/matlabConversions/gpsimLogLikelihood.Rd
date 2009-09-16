\name{gpsimLogLikelihood}
\Rdversion{1.0}
\alias{gpsimLogLikelihood}
\title{Compute the log likelihood of a GPSIM model.}
\description{
  computes the log likelihood of the given Gaussian process for use in a single input motif protein network.
}
\usage{
ll <- gpsimLogLikelihood(model)
}
\arguments{
  \item{model}{the model for which the log likelihood is computed.}
}
\value{
  \item{ll}{the log likelihood of the data set. }
}
\seealso{
\code{\link{gpsimCreate}, \link{gpsimLogLikeGradient}, \link{gpsimObjective}}.
}
\examples{
## missing
}
\keyword{model}
