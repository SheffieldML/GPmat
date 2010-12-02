\name{noiseOut}
\Rdversion{1.0}
\alias{noiseOut}
\title{Give the output of the noise model given the mean and variance.}
\description{
  computes the ouptut for the given noise given input mean and variances.
}
\usage{
y <- noiseOut(noise, mu, varSigma)
}
\arguments{
  \item{noise}{the noise structure for which the output is computed.}
  \item{mu}{the input mean values.}
  \item{varSigma}{the input variance values.}
}
\value{
  \item{y}{the output from the noise model.}
}
\seealso{
\code{\link{noiseParamInit}, \link{noiseCreate}}.
}
\examples{
## missing
}
\keyword{model}
