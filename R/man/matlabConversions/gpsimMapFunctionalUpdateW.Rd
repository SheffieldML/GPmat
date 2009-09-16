\name{gpsimMapFunctionalUpdateW}
\Rdversion{1.0}
\alias{gpsimMapFunctionalUpdateW}
\title{Update the data component of the Hessian.}
\description{
  updates the portion of the Hessian associated with the likelihood of the data.
}
\usage{
model <- gpsimMapFunctionalUpdateW(model)
}
\arguments{
  \item{model}{the model for which the data component of the Hessian is to be updated.}
}
\value{
  \item{model}{the model with the data component of the Hessian up to date.}
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapFunctionalGradient}, \link{gpsimMapUpdatePosteriorCovariance}}.
}
\examples{
## missing
}
\keyword{model}
