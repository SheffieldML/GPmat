\name{translateKernExtractParam}
\Rdversion{1.0}
\alias{translateKernExtractParam}
\title{Extract parameters from the TRANSLATE kernel structure.}
\description{
  Extract parameters from the TRANSLATE kernel structure.
}
\usage{
param <- translateKernExtractParam(kern)
values <- translateKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\details{
  \code{param <- translateKernExtractParam(kern)}
  extracts parameters from the input space translation kernel structure into a vector of parameters for optimisation.

  \code{param, names <- translateKernExtractParam(kern)}
  extracts parameters and parameter names from the input space translation kernel structure.

}
\value{
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel structure, the parameters will be transformed before optimisation (for example positive only parameters could be logged before being returned).}
  \item{names}{cell array of strings containing names for each parameter.}
}
\seealso{
\code{\link{kernExtractParam}, \link{cmpndKernExtractParam}, \link{scg}, \link{conjgrad}}.
}
\examples{
## missing
}
\keyword{model}
