\name{multiKernDiagCompute}
\Rdversion{1.0}
\alias{multiKernDiagCompute}
\title{Compute diagonal of MULTI kernel.}
\description{
  computes the diagonal of the kernel matrix for the multiple output block kernel given a design matrix of inputs.
}
\usage{
k <- multiKernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{multiKernCompute}}.
}
\examples{
## missing
}
\keyword{model}
