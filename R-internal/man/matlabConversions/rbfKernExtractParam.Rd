\name{rbfKernExtractParam}
\Rdversion{1.0}
\alias{rbfKernExtractParam}
\title{Extract parameters from the RBF kernel structure.}
\description{
  Extract parameters from the RBF kernel structure.
}
\usage{
param <- rbfKernExtractParam(kern)
values <- rbfKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\details{
  \code{param <- rbfKernExtractParam(kern)}
  Extract parameters from the radial basis function kernel structure into a vector of parameters for optimisation.

  \code{param, names <- rbfKernExtractParam(kern)}
  Extract parameters and parameter names from the radial basis function kernel structure.

}
\value{
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel matrix, the parameters will be transformed before optimisation (for example positive only parameters could be logged before being returned).}
  \item{names}{cell array of strings giving names to the parameters.}
}
\seealso{
\code{\link{}}.
}
\examples{
## missing
}
\keyword{model}
