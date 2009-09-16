\name{cmpndKernDiagGradX}
\Rdversion{1.0}
\alias{cmpndKernDiagGradX}
\title{Gradient of CMPND kernel's diagonal with respect to X.}
\description{
  computes the gradient of the diagonal of the compound kernel matrix with respect to the elements of the design matrix given in X.
}
\usage{
gX <- cmpndKernDiagGradX(kern, X)
}
\arguments{
  \item{kern}{the kernel structure for which gradients are being computed.}
  \item{X}{the input data in the form of a design matrix.}
}
\value{
  \item{gX}{the gradients of the diagonal with respect to each element of X. The returned matrix has the same dimensions as X.}
}
\seealso{
\code{\link{cmpndKernParamInit}, \link{kernDiagGradX}, \link{cmpndkernGradX}}.
}
\examples{
## missing
}
\keyword{model}
