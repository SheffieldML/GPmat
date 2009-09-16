\name{multiKernExpandParam}
\Rdversion{1.0}
\alias{multiKernExpandParam}
\title{Create kernel structure from MULTI kernel's parameters.}
\description{
  returns a multiple output block kernel structure filled with the parameters in the given vector. This is used as a helper function to enable parameters to be optimised in, for example, the NETLAB optimisation functions.
}
\usage{
kern <- multiKernExpandParam(kern, param)
}
\arguments{
  \item{kern}{the kernel structure in which the parameters are to be placed.}
  \item{param}{vector of parameters which are to be placed in the kernel structure.}
}
\value{
  \item{kern}{kernel structure with the given parameters in the relevant locations.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernExtractParam}, \link{kernExpandParam}}.
}
\examples{
## missing
}
\keyword{model}
