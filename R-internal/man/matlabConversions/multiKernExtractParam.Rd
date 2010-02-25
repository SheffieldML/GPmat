\name{multiKernExtractParam}
\Rdversion{1.0}
\alias{multiKernExtractParam}
\title{Extract parameters from the MULTI kernel structure.}
\description{
  Extract parameters from the MULTI kernel structure.
}
\usage{
param <- multiKernExtractParam(kern)
values <- multiKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\details{
  \code{param <- multiKernExtractParam(kern)}
  Extract parameters from the multiple output block kernel structure into a vector of parameters for optimisation.

  \code{param, names <- multiKernExtractParam(kern)}
  Extract parameters and parameter names from the multiple output block kernel structure.

}
\value{
  \item{param}{vector of parameters extracted from the kernel. The vector of transforms is assumed to be empty here, any transormation of parameters is assumed to be done in the component kernels.}
  \item{names}{cell array of strings containing parameter names.}
}
\seealso{
\code{\link{}}.
}
\examples{
## missing
}
\keyword{model}
