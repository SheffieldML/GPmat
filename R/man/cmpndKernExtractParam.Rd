\name{cmpndKernExtractParam}
\Rdversion{1.0}
\alias{cmpndKernExtractParam}
\title{Extract parameters from the CMPND kernel structure.}
\description{
  Extract parameters from the compound kernel matrix into a vector of parameters for optimisation.
}
\usage{
param <- cmpndKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\value{
  \item{param}{vector of parameters extracted from the kernel. The vector of 'transforms' is assumed to be empty here. Any transformations of parameters should be done in component kernels.}
}
\seealso{
\code{\link{}}.
}
\examples{
## missing
}
\keyword{model}
