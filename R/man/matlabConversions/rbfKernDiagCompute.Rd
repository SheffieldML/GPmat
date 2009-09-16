\name{rbfKernDiagCompute}
\Rdversion{1.0}
\alias{rbfKernDiagCompute}
\title{Compute diagonal of RBF kernel.}
\description{
  computes the diagonal of the kernel matrix for the radial basis function kernel given a design matrix of inputs.
}
\usage{
k <- rbfKernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{rbfKernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{rbfKernCompute}}.
}
\examples{
## missing
}
\keyword{model}
