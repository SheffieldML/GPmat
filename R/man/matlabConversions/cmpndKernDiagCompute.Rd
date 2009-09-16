\name{cmpndKernDiagCompute}
\Rdversion{1.0}
\alias{cmpndKernDiagCompute}
\title{Compute diagonal of CMPND kernel.}
\description{
  computes the diagonal of the kernel matrix for the compound kernel given a design matrix of inputs.
}
\usage{
k <- cmpndKernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{cmpndKernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{cmpndKernCompute}}.
}
\examples{
## missing
}
\keyword{model}
