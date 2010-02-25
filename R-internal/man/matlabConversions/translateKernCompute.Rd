\name{translateKernCompute}
\Rdversion{1.0}
\alias{translateKernCompute}
\title{Compute the TRANSLATE kernel given the parameters and X.}
\description{
  Compute the TRANSLATE kernel given the parameters and X.
}
\usage{
k <- translateKernCompute(kern, x, x2)
k <- translateKernCompute(kern, x')
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{the input matrix associated with the rows of the kernel.}
  \item{x2}{the input matrix associated with the columns of the kernel.}
  \item{x'}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k <- translateKernCompute(kern, x, x2)}
  computes the kernel parameters for the input space translation kernel given inputs associated with rows and columns.

  \code{k <- translateKernCompute(kern, x')}
  computes the kernel matrix for the input space translation kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{translateKernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{cmpndKernCompute}, \link{translateKernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
