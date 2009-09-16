\name{cgpsimLogLikelihood}
\Rdversion{1.0}
\alias{cgpsimLogLikelihood}
\title{Compound GPSIM model's log likelihood.}
\description{
  this is just a wrapper for learning three GP models simultaneously.
}
\usage{
ll <- cgpsimLogLikelihood(model)
}
\arguments{
  \item{model}{input model for which gradients of the likelihood will be computed.}
}
\value{
  \item{ll}{the log likelihood.}
}
\seealso{
\code{\link{gpsimLogLikelihood}}.
}
\examples{
## missing
}
\keyword{model}
