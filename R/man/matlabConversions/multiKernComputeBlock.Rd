\name{multiKernComputeBlock}
\Rdversion{1.0}
\alias{multiKernComputeBlock}
\title{}
\description{
  
}
\usage{
K <- multiKernComputeBlock(kern, X1, X2, i, j)
K <- multiKernComputeBlock(kern, X, i, j)
}
\arguments{
  \item{kern}{the structure containing the kernel.}
  \item{X1}{first set of kernel inputs.}
  \item{X2}{second set of kernel inputs.}
  \item{i}{the row of the block of the kernel to be computed.}
  \item{j}{the column of the block of the kernel to be computed.}
  \item{X}{first set of kernel inputs.}
}
\details{
  \code{K <- multiKernComputeBlock(kern, X1, X2, i, j)}
  computes a block of a multi-output kernel given two matrices of input.

  \code{K <- multiKernComputeBlock(kern, X, i, j)}
  computes a block of a multi-output kernel given a single matrix of input.

}
\value{
  \item{K}{the kernel matrix for the given inputs.}
}
\seealso{
\code{\link{multiKernCreate}, \link{multiKernCompute}, \link{multiKernGradientBlock}}.
}
\examples{
## missing
}
\keyword{model}
