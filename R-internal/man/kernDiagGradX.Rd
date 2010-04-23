\name{kernDiagGradX}
\Rdversion{1.0}
\alias{kernGradX}
\alias{kernDiagGradX}
\alias{cmpndKernGradX}
\alias{cmpndKernDiagGradX}
\alias{mlpKernDiagGradX}
\alias{mlpKernGradX}
\title{Compute the gradient of the  kernel wrt X.}
\description{
  computes the gradient of the (diagonal of the) kernel matrix with respect to the elements of the design matrix given in X.
}
\usage{
gX <- kernDiagGradX(kern, X)
gX2 <- kernGradX(kern, X)
gX2 <- kernGradX(kern, X1, X2)
}
\arguments{
  \item{kern}{the kernel structure for which gradients are being computed.}
  \item{X}{the input data in the form of a design matrix.}
  \item{X1}{row locations against which gradients are being computed.}
  \item{X2}{column locations against which gradients are being computed.}
}
\value{
  \item{gX}{the gradients of the diagonal with respect to each element of X. The returned matrix has the same dimensions as X.}
  \item{gX2}{the returned gradients. The gradients are returned in a matrix which is numData x numInputs x numData. Where numData is the number of data points and numInputs is the number of input dimensions in X.}
}
\seealso{
\code{\link{kernGradient}}
}
\examples{
kern <- kernCreate(1, 'mlp')
g <- kernDiagGradX(kern, as.matrix(3:8))
}
\keyword{model}
