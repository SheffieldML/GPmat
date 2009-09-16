\name{kernCompute}
\Rdversion{1.0}
\alias{kernCompute}
\title{Compute the kernel given the parameters and X.}
\description{
  Compute the kernel given the parameters and X.
}
\usage{
K <- kernCompute(kern, X)
K <- kernCompute(kern, X', X2)
}
\arguments{
  \item{kern}{kernel structure to be computed.}
  \item{X}{input data matrix (rows are data points) to the kernel computation.}
  \item{X'}{first input matrix to the kernel computation (forms the rows of the kernel).}
  \item{X2}{second input matrix to the kernel computation (forms the columns of the kernel).}
}
\details{
  \code{K <- kernCompute(kern, X)}
  computes a kernel matrix for the given kernel type given an input data matrix.

  \code{K <- kernCompute(kern, X', X2)}
  computes a kernel matrix for the given kernel type given two input data matrices, one for the rows and one for the columns.

}
\value{
  \item{K}{computed elements of the kernel structure.}
}
\seealso{
\code{\link{kernCreate}, \link{kernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
