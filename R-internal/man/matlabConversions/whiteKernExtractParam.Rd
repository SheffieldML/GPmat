\name{whiteKernExtractParam}
\Rdversion{1.0}
\alias{whiteKernExtractParam}
\title{Extract parameters from the WHITE kernel structure.}
\description{
  Extract parameters from the WHITE kernel structure.
}
\usage{
param <- whiteKernExtractParam(kern)
values <- whiteKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\details{
  \code{param <- whiteKernExtractParam(kern)}
  Extract parameters from the white noise kernel structure into a vector of parameters for optimisation.

  \code{param, names <- whiteKernExtractParam(kern)}
  Extract parameters and parameter names from the white noise kernel structure.

}
\value{
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel matrix, the parameters will be transformed before optimisation (for example positive only parameters could be logged before being returned).}
  \item{names}{cell array of strings giving paramter names.}
}
\seealso{
\code{\link{}}.
}
\examples{
## missing
}
\keyword{model}
