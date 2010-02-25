\name{disimKernExpandParam}
\Rdversion{1.0}
\alias{disimKernExpandParam}
\title{Create kernel structure from DISIM kernel's parameters.}
\description{
  returns a single input motif kernel structure filled with the parameters in the given vector. This is used as a helper function to enable parameters to be optimised in, for example, the NETLAB optimisation functions.
}
\usage{
kern <- disimKernExpandParam(kern, param)
}
\arguments{
  \item{kern}{the kernel structure in which the parameters are to be placed.}
  \item{param}{vector of parameters which are to be placed in the kernel structure.}
}
\value{
  \item{kern}{kernel structure with the given parameters in the relevant locations.}
}
\seealso{
\code{\link{disimKernParamInit}, \link{disimKernExtractParam}, \link{kernExpandParam}}.
}
\examples{
## missing
}
\keyword{model}
