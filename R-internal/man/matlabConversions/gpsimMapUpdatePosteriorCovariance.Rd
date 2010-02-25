\name{gpsimMapUpdatePosteriorCovariance}
\Rdversion{1.0}
\alias{gpsimMapUpdatePosteriorCovariance}
\title{update the posterior covariance of f.}
\description{
  updates the stored representation of the posterior covariance of f in a GPSIMMAP model. Relies on the representation of the kernel and the data Hessian being up to date. These can be forced to be up to date by using gpsimMapUpdateKernels and gpsimMapFunctionalUpdateW.
}
\usage{
model <- gpsimMapUpdatePosteriorCovariance(model)
}
\arguments{
  \item{model}{the model for which the representation is to be updated.}
}
\value{
  \item{model}{the model with the updated representation.}
}
\seealso{
\code{\link{gpsimMapUpdateKernels}, \link{gpsimMapFunctionalUpdateW}}.
}
\examples{
## missing
}
\keyword{model}
