\name{simKernDiagCompute}
\Rdversion{1.0}
\alias{simKernDiagCompute}
\title{Compute diagonal of SIM kernel.}
\description{
  computes the diagonal of the kernel matrix for the single input motif kernel given a design matrix of inputs.
}
\usage{
k <- simKernDiagCompute(kern, t)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{t}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{simKernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{simKernCompute}}.
}
\examples{
## missing
}
\keyword{model}
