\name{gpsimMapObjective}
\Rdversion{1.0}
\alias{gpsimMapObjective}
\title{Compute the objective function  of the GPSIMMAP model.}
\description{
  computes the objective function of the model given Gaussian process for use in a single input motif protein network.
}
\usage{
f <- gpsimMapObjective(model)
}
\arguments{
  \item{model}{the model for which the objective function is computed.}
}
\value{
  \item{f}{the objective function of the model. }
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapLogLikelihood}, \link{gpsimMapGradient}, \link{gpsimMapFunctionalLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
