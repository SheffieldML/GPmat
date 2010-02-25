\name{rbfKernGradient}
\Rdversion{1.0}
\alias{rbfKernGradient}
\title{Gradient of RBF kernel's parameters.}
\description{
  Gradient of RBF kernel's parameters.
}
\usage{
g <- rbfKernGradient(kern, x, partial)
g <- rbfKernGradient(kern, x1, x2, partial')
}
\arguments{
  \item{kern}{the kernel structure for which the gradients are being computed.}
  \item{x}{the input locations for which the gradients are being computed.}
  \item{partial}{matrix of partial derivatives of the function of interest with respect to the kernel matrix. The argument takes the form of a square matrix of dimension  numData, where numData is the number of rows in X.}
  \item{x1}{the input locations associated with the rows of the kernel matrix.}
  \item{x2}{the input locations associated with the columns of the kernel matrix.}
  \item{partial'}{matrix of partial derivatives of the function of interest with respect to the kernel matrix. The matrix should have the same number of rows as X1 and the same number of columns as X2 has rows.}
}
\details{
  \code{g <- rbfKernGradient(kern, x, partial)}
  computes the gradient of functions with respect to the radial basis function kernel's parameters. As well as the kernel structure and the input positions, the user provides a matrix PARTIAL which gives the partial derivatives of the function with respect to the relevant elements of the kernel matrix.

  \code{g <- rbfKernGradient(kern, x1, x2, partial')}
  computes the derivatives as above, but input locations are now provided in two matrices associated with rows and columns of the kernel matrix.

}
\value{
  \item{g}{gradients of the function of interest with respect to the kernel parameters. The ordering of the vector should match that provided by the function kernExtractParam.}
}
\seealso{
\code{\link{}}.
}
\examples{
## missing
}
\keyword{model}
