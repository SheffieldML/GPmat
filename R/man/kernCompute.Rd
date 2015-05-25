\name{kernCompute}
\Rdversion{1.0}
\alias{kernCompute}
\alias{kernDiagCompute}
\alias{cmpndKernCompute}
\alias{cmpndKernDiagCompute}
\alias{disimKernCompute}
\alias{disimKernDiagCompute}
\alias{mlpKernCompute}
\alias{multiKernCompute}
\alias{multiKernDiagCompute}
\alias{rbfKernCompute}
\alias{rbfKernDiagCompute}
\alias{simKernCompute}
\alias{simKernDiagCompute}
\alias{translateKernCompute}
\alias{translateKernDiagCompute}
\alias{whiteKernCompute}
\alias{whiteKernDiagCompute}
\alias{disimXdisimKernCompute}
\alias{disimXrbfKernCompute}
\alias{disimXsimKernCompute}
\alias{simXrbfKernCompute}
\alias{simXsimKernCompute}
\alias{whiteXwhiteKernCompute}
\title{Compute the kernel given the parameters and X.}
\description{
  Compute the kernel given the parameters and X.
}
\usage{
  kernCompute(kern, x, x2)
  kernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{kernel structure to be computed.}
  \item{x}{first or only input data matrix (rows are data points) to the kernel computation.}
  \item{x2}{(optional) second input matrix to the kernel computation (forms the columns of the kernel).}
}
\details{
  \code{K <- kernCompute(kern, X)}
  computes a kernel matrix for the given kernel type given an input data matrix.

  \code{K <- kernCompute(kern, X1, X2)}
  computes a kernel matrix for the given kernel type given two input data matrices, one for the rows and one for the columns.
  
  \code{K <- kernDiagCompute(kern, X)}
  computes the diagonal of a kernel matrix for the given kernel.
  
  \code{K <- *X*kernCompute(kern1, kern2, X)}
  \code{K <- *X*kernCompute(kern1, kern2, X1, X2)}
  same as above, but for cross combinations of two kernels, \code{kern1} and \code{kern2}.
}
\value{
  \item{K}{computed elements of the kernel structure.}
  \item{Kd}{vector containing computed diagonal elements of the kernel structure.}
}
\seealso{
\code{\link{kernCreate}}
}
\examples{
kern <- kernCreate(1, 'rbf')
K <- kernCompute(kern, as.matrix(3:8))
}
\keyword{model}
