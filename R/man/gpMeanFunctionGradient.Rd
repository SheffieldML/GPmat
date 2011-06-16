\name{gpMeanFunctionGradient}
\Rdversion{1.0}
\alias{gpMeanFunctionGradient}
\title{Compute the log likelihood gradient wrt the scales.}
\description{
  computes the gradient of the log likelihood with respect to the scales. In the future the gradients with respect to the biases may also be included.
}
\usage{
  gpMeanFunctionGradient(model)
}
\arguments{
  \item{model}{the model for which the gradients are to be computed.}
}
\value{
  \item{g}{the gradients of the likelihood with respect to the mean function's parameters.}
}
\seealso{
\code{\link{gpCreate}, \link{gpScaleBiasGradient}, \link{gpLogLikeGradients}, \link{gpLogLikelihood}}.
}
\examples{
## missing
}
\keyword{model}
