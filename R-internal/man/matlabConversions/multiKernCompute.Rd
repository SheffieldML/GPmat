\name{multiKernCompute}
\Rdversion{1.0}
\alias{multiKernCompute}
\title{Compute the MULTI kernel given the parameters and X.}
\description{
  Compute the MULTI kernel given the parameters and X.
}
\usage{
k <- multiKernCompute(kern, x, x2)
k <- multiKernCompute(kern, x')
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{the input matrix associated with the rows of the kernel.}
  \item{x2}{the input matrix associated with the columns of the kernel.}
  \item{x'}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k <- multiKernCompute(kern, x, x2)}
  computes the kernel parameters for the multiple output block kernel given inputs associated with rows and columns.

  \code{k <- multiKernCompute(kern, x')}
  computes the kernel matrix for the multiple output block kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{multiKernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
