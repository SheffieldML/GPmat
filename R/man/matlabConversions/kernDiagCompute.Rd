\name{kernDiagCompute}
\Rdversion{1.0}
\alias{kernDiagCompute}
\title{Compute the kernel given the parameters and X.}
\description{
  computes the diagonal of a kernel matrix for the given kernel.
}
\usage{
K <- kernDiagCompute(kern, X)
}
\arguments{
  \item{kern}{kernel structure to be computed.}
  \item{X}{input data matrix (rows are data points) to the kernel computation.}
}
\value{
  \item{K}{vector containing computed diagonal elements of the kernel structure.}
}
\seealso{
\code{\link{kernCompute}}.
}
\examples{
## missing
}
\keyword{model}
