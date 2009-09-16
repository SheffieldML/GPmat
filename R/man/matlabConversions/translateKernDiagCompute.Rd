\name{translateKernDiagCompute}
\Rdversion{1.0}
\alias{translateKernDiagCompute}
\title{Compute diagonal of TRANSLATE kernel.}
\description{
  computes the diagonal of the kernel matrix for the input space translation kernel given a design matrix of inputs.
}
\usage{
k <- translateKernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{translateKernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{cmpndKernDiagCompute}, \link{translateKernCompute}}.
}
\examples{
## missing
}
\keyword{model}
