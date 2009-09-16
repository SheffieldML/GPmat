\name{multiKernGradientBlock}
\Rdversion{1.0}
\alias{multiKernGradientBlock}
\title{}
\description{
  
}
\usage{
values <- multiKernGradientBlock(kern, X1, X2, covGrad, i, j)
values <- multiKernGradientBlock(kern, X, covGrad, i, j)
}
\arguments{
  \item{kern}{the structure containing the kernel.}
  \item{X1}{first set of kernel inputs.}
  \item{X2}{second set of kernel inputs.}
  \item{covGrad}{Gradient of the objective function with respect to the relevant portion of the kernel matrix.}
  \item{i}{the row of the block of the kernel to be computed.}
  \item{j}{the column of the block of the kernel to be computed.}
  \item{X}{first set of kernel inputs.}
}
\details{
  \code{g1, g2 <- multiKernGradientBlock(kern, X1, X2, covGrad, i, j)}
  computes the gradient with respect to parameters for a block of a multi-output kernel given two matrices of input.

  \code{g1, g2 <- multiKernGradientBlock(kern, X, covGrad, i, j)}
  compute a block of a multi-output kernel given a single matrix of input.

}
\value{
  \item{g1}{the gradient of the kernel parameters from the first kernel in the order provided by the relevant kernExtractParam commands.}
  \item{g2}{the gradient of the kernel parameters from the second kernel in the order provided by the relevant kernExtractParam commands.}
}
\seealso{
\code{\link{multiKernCreate}, \link{multiKernGradient}, \link{multiKernComputeBlock}}.
}
\examples{
## missing
}
\keyword{model}
