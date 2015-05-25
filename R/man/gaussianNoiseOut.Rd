\name{gaussianNoiseOut}
\Rdversion{1.0}
\alias{gaussianNoiseOut}
\title{Compute the output of the GAUSSIAN noise given the input mean and variance.}
\description{
  computes the ouptut for the Gaussian noise given input mean and variances.
}
\usage{
gaussianNoiseOut(noise, mu, varSigma)
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
\code{\link{gaussianNoiseParamInit}, \link{noiseOut}, \link{noiseCreate}}.
}
\examples{
## missing
}
\keyword{model}
