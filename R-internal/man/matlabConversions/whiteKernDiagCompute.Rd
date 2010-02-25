\name{whiteKernDiagCompute}
\Rdversion{1.0}
\alias{whiteKernDiagCompute}
\title{Compute diagonal of WHITE kernel.}
\description{
  computes the diagonal of the kernel matrix for the white noise kernel given a design matrix of inputs.
}
\usage{
k <- whiteKernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{whiteKernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{whiteKernCompute}}.
}
\examples{
## missing
}
\keyword{model}
