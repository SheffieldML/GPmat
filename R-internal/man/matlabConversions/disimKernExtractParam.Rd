\name{disimKernExtractParam}
\Rdversion{1.0}
\alias{disimKernExtractParam}
\title{Extract parameters from the DISIM kernel structure.}
\description{
  Extract parameters from the DISIM kernel structure.
}
\usage{
param <- disimKernExtractParam(kern)
values <- disimKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\details{
  \code{param <- disimKernExtractParam(kern)}
  Extract parameters from the single input motif kernel structure into a vector of parameters for optimisation.

  \code{param, names <- disimKernExtractParam(kern)}
  Extract parameters and their names from the single input motif kernel structure.

}
\value{
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel matrix, the parameters will be transformed before optimisation (for example positive only parameters could be logged before being returned).}
  \item{names}{cell array of strings containing parameter names.}
}
\seealso{
\code{\link{}}.
}
\examples{
## missing
}
\keyword{model}
