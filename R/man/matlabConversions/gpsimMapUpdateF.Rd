\name{gpsimMapUpdateF}
\Rdversion{1.0}
\alias{gpsimMapUpdateF}
\title{Update posterior mean of f.}
\description{
  updates the stored version of the posterior mean for f in the GPSIMMAP model.
}
\usage{
values <- gpsimMapUpdateF(model, options)
}
\arguments{
  \item{model}{the model for which f is to be updated.}
  \item{options}{options vector.}
}
\value{
  \item{model}{the model with f updated.}
  \item{ll}{the log likelihood after update.}
}
\seealso{
\code{\link{gpsimMapFunctionalLogLikeGradient}, \link{gpsimMapFunctionalLogLikeHessian}, \link{gpsimMapCreate}}.
}
\examples{
## missing
}
\keyword{model}
