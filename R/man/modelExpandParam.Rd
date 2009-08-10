\name{modelExpandParam}
\Rdversion{1.0}
\alias{modelExpandParam}
\title{Update a model structure with parameters.}
\description{
  returns a model structure filled with the parameters in the given vector. This is used as a helper function to enable parameters to be optimised in, for example, the NETLAB optimisation functions.
}
\usage{
model <- modelExpandParam(model, param)
}
\arguments{
  \item{model}{the model structure in which the parameters are to be placed.}
  \item{param}{vector of parameters which are to be placed in the model structure.}
}
\value{
  \item{model}{model structure with the given parameters in the relevant locations.}
}
\seealso{
\code{\link{modelExtractParam}, \link{scg}, \link{conjgrad}}.
}
\examples{
## missing
}
\keyword{model}
