\name{mlpKernExtractParam}
\Rdversion{1.0}
\alias{mlpKernExtractParam}
\title{Extract parameters from the MLP kernel structure.}
\description{
  Extract parameters from the MLP kernel structure.
}
\usage{
param <- mlpKernExtractParam(kern)
values <- mlpKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\details{
  \code{param <- mlpKernExtractParam(kern)}
  Extract parameters from the multi-layer perceptron kernel matrix into a vector of parameters for optimisation.

  \code{param, names <- mlpKernExtractParam(kern)}
  Extract parameters and parameter names from the multi-layer perceptron kernel matrix into a vector of parameters for optimisation.

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
