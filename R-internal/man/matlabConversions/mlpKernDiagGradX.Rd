\name{mlpKernDiagGradX}
\Rdversion{1.0}
\alias{mlpKernDiagGradX}
\title{Gradient of MLP kernel's diagonal with respect to X.}
\description{
  computes the gradient of the diagonal of the multi-layer perceptron kernel matrix with respect to the elements of the design matrix given in X.
}
\usage{
gX <- mlpKernDiagGradX(kern, X)
}
\arguments{
  \item{kern}{the kernel structure for which gradients are being computed.}
  \item{X}{the input data in the form of a design matrix.}
}
\value{
  \item{gX}{the gradients of the diagonal with respect to each element of X. The returned matrix has the same dimensions as X.}
}
\seealso{
\code{\link{mlpKernParamInit}, \link{kernDiagGradX}, \link{mlpkernGradX}}.
}
\examples{
## missing
}
\keyword{model}
