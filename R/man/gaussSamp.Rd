\name{gaussSamp}
\Rdversion{1.0}
\alias{gaussSamp}
\title{Sample from a Gaussian with a given covariance.}
\description{
  samples a given number of samples from a Gaussian with a given covariance matrix.
}
\usage{
gaussSamp(mu=matrix(0,nrow=dim(Sigma)[1]), Sigma, numSamps)
}
\arguments{
  \item{mu}{the mean of the Gaussian to sample from.}
  \item{Sigma}{the covariance of the Gaussian to sample from.}
  \item{numSamps}{the number of samples to take from the Gaussian.}
}
\value{
  \item{y}{the samples from the Gaussian}
}
\seealso{
\code{\link{rnorm}, \link{eigen}}.
}
\examples{
## missing
}
\keyword{model}
