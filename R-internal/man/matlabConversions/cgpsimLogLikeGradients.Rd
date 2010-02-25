\name{cgpsimLogLikeGradients}
\Rdversion{1.0}
\alias{cgpsimLogLikeGradients}
\title{Compound GPSIM model's gradients.}
\description{
  this is just a wrapper for learning three GP models simultaneously.
}
\usage{
g <- cgpsimLogLikeGradients(model)
}
\arguments{
  \item{model}{input model for which gradients of the likelihood will be computed.}
}
\value{
  \item{g}{gradients of the log likelihood with repsect to the parameters.}
}
\seealso{
\code{\link{gpsimLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
