\name{disimKernDiagCompute}
\Rdversion{1.0}
\alias{disimKernDiagCompute}
\title{Compute diagonal of DISIM kernel.}
\description{
  computes the diagonal of the kernel matrix for the driven input single input motif kernel given a design matrix of inputs.
}
\usage{
k <- disimKernDiagCompute(kern, t)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{t}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{disimKernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{disimKernCompute}}.
}
\examples{
## missing
}
\keyword{model}
