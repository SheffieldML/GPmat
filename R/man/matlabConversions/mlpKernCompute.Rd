\name{mlpKernCompute}
\Rdversion{1.0}
\alias{mlpKernCompute}
\title{Compute the MLP kernel given the parameters and X.}
\description{
  Compute the MLP kernel given the parameters and X.
}
\usage{
k <- mlpKernCompute(kern, x, x2)
k <- mlpKernCompute(kern, x')
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{the input matrix associated with the rows of the kernel.}
  \item{x2}{the input matrix associated with the columns of the kernel.}
  \item{x'}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k <- mlpKernCompute(kern, x, x2)}
  computes the kernel parameters for the multi-layer perceptron kernel given inputs associated with rows and columns.

  \code{k <- mlpKernCompute(kern, x')}
  computes the kernel matrix for the multi-layer perceptron kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{mlpKernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{mlpKernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
