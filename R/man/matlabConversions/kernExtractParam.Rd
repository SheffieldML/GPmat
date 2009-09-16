\name{kernExtractParam}
\Rdversion{1.0}
\alias{kernExtractParam}
\title{Extract parameters from kernel structure.}
\description{
  Extract parameters from the kernel into a vector of parameters for optimisation.
}
\usage{
param <- kernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\value{
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel, the parameters will be transformed before optimisation (for example positive only parameters could be logged before being returned).}
}
\seealso{
\code{\link{kernExpandParam}, \link{scg}, \link{conjgrad}}.
}
\examples{
## missing
}
\keyword{model}
