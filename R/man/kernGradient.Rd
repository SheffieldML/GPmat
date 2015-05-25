\name{kernGradient}
\Rdversion{1.0}
\alias{kernGradient}
\alias{cmpndKernGradient}
\alias{disimKernGradient}
\alias{disimXdisimKernGradient}
\alias{disimXrbfKernGradient}
\alias{disimXsimKernGradient}
\alias{mlpKernGradient}
\alias{multiKernGradient}
\alias{rbfKernGradient}
\alias{simKernGradient}
\alias{simXrbfKernGradient}
\alias{simXsimKernGradient}
\alias{translateKernGradient}
\alias{whiteKernGradient}
\alias{whiteXwhiteKernGradient}
\title{Compute the gradient wrt the kernel parameters.}
\description{
  Compute the gradient wrt the kernel parameters.
}
\usage{
  kernGradient(kern, x, ...)
}
\arguments{
  \item{kern}{the kernel structure for which the gradients are being computed.}
  \item{x}{the input locations for which the gradients are being computed,
    specifically those associated with the rows of the kernel matrix
    if there are two arguments of input locations.}
  \item{...}{optional arguments including potentially:
    the input locations associated with the columns of the kernel matrix;
    matrix of partial derivatives of the function of interest with
    respect to the kernel matrix. With single input, the argument
    takes the form of a square matrix of dimension numData, where numData
    is the number of rows in x, with two input arguments
    the matrix should have the same number of rows as the first and the same
    number of columns as the second has rows.
    }
}
\details{
  \code{g <- kernGradient(kern, x, partial)}
  \code{g <- *kernGradient(kern, x, partial)}
  computes the gradient of functions with respect to the kernel parameters. As well as the kernel structure and the input positions, the user provides a matrix PARTIAL which gives the partial derivatives of the function with respect to the relevant elements of the kernel matrix.

  \code{g <- kernGradient(kern, x1, x2, partial_)}
  \code{g <- *kernGradient(kern, x1, x2, partial_)}
  computes the derivatives as above, but input locations are now provided in two matrices associated with rows and columns of the kernel matrix.

  \code{g <- *X*kernGradient(kern1, kern2, x, partial)}
  \code{g <- *X*kernGradient(kern1, kern2, x1, x2, partial_)}
  same as above, but for cross combinations of two kernels, \code{kern1} and \code{kern2}.
}
\value{
  \item{g}{gradients of the function of interest with respect to the kernel parameters. The ordering of the vector should match that provided by the function kernExtractParam.}
}
\seealso{
\code{\link{kernCompute}}, \code{\link{kernExtractParam}}.
}
\examples{
kern <- kernCreate(1, 'rbf')
g <- kernGradient(kern, as.matrix(c(1, 4)), array(1, c(2, 2)))
}
\keyword{model}
