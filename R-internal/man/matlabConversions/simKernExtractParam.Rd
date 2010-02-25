\name{simKernExtractParam}
\Rdversion{1.0}
\alias{simKernExtractParam}
\title{Extract parameters from the SIM kernel structure.}
\description{
  Extract parameters from the SIM kernel structure.
}
\usage{
param <- simKernExtractParam(kern)
values <- simKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\details{
  \code{param <- simKernExtractParam(kern)}
  Extract parameters from the single input motif kernel structure into a vector of parameters for optimisation.

  \code{param, names <- simKernExtractParam(kern)}
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
