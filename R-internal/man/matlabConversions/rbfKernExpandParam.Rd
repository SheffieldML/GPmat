\name{rbfKernExpandParam}
\Rdversion{1.0}
\alias{rbfKernExpandParam}
\title{Create kernel structure from RBF kernel's parameters.}
\description{
  returns a radial basis function kernel structure filled with the parameters in the given vector. This is used as a helper function to enable parameters to be optimised in, for example, the NETLAB optimisation functions.
}
\usage{
kern <- rbfKernExpandParam(kern, param)
}
\arguments{
  \item{kern}{the kernel structure in which the parameters are to be placed.}
  \item{param}{vector of parameters which are to be placed in the kernel structure.}
}
\value{
  \item{kern}{kernel structure with the given parameters in the relevant locations.}
}
\seealso{
\code{\link{rbfKernParamInit}, \link{rbfKernExtractParam}, \link{kernExpandParam}}.
}
\examples{
## missing
}
\keyword{model}
