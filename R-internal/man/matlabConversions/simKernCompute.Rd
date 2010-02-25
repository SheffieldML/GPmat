\name{simKernCompute}
\Rdversion{1.0}
\alias{simKernCompute}
\title{Compute the SIM kernel given the parameters and X.}
\description{
  Compute the SIM kernel given the parameters and X.
}
\usage{
values <- simKernCompute(kern, t1, t2)
values <- simKernCompute(kern, t)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{t1}{the input matrix associated with the rows of the kernel.}
  \item{t2}{the input matrix associated with the columns of the kernel.}
  \item{t}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k, sk <- simKernCompute(kern, t1, t2)}
  computes the kernel parameters for the single input motif kernel given inputs associated with rows and columns.

  \code{k, sk <- simKernCompute(kern, t)}
  computes the kernel matrix for the single input motif kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
  \item{sk}{unscaled kernel matrix (i.e. only 0.5 times the sum of h's part).}
}
\seealso{
\code{\link{simKernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{simKernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
