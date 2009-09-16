\name{disimKernCompute}
\Rdversion{1.0}
\alias{disimKernCompute}
\title{Compute the DISIM kernel given the parameters and X.}
\description{
  Compute the DISIM kernel given the parameters and X.
}
\usage{
k <- disimKernCompute(kern, t1, t2)
k <- disimKernCompute(kern, t)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{t1}{the input matrix associated with the rows of the kernel.}
  \item{t2}{the input matrix associated with the columns of the kernel.}
  \item{t}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k <- disimKernCompute(kern, t1, t2)}
  computes the kernel parameters for the single input motif kernel given inputs associated with rows and columns.

  \code{k <- disimKernCompute(kern, t)}
  computes the kernel matrix for the single input motif kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{disimKernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{disimKernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
