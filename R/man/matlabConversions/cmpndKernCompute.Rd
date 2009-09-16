\name{cmpndKernCompute}
\Rdversion{1.0}
\alias{cmpndKernCompute}
\title{Compute the CMPND kernel given the parameters and X.}
\description{
  Compute the CMPND kernel given the parameters and X.
}
\usage{
k <- cmpndKernCompute(kern, x, x2)
k <- cmpndKernCompute(kern, x')
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{the input matrix associated with the rows of the kernel.}
  \item{x2}{the input matrix associated with the columns of the kernel.}
  \item{x'}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k <- cmpndKernCompute(kern, x, x2)}
  computes the kernel parameters for the compound kernel given inputs associated with rows and columns.

  \code{k <- cmpndKernCompute(kern, x')}
  computes the kernel matrix for the compound kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{cmpndKernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{cmpndKernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
