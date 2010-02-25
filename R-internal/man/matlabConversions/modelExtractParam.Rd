\name{modelExtractParam}
\Rdversion{1.0}
\alias{modelExtractParam}
\title{Extract the parameters of a model.}
\description{
  Extract parameters from the model into a vector of parameters for optimisation.
}
\usage{
param <- modelExtractParam(model)
}
\arguments{
  \item{model}{the model structure containing the parameters to be extracted.}
}
\value{
  \item{param}{vector of parameters extracted from the model.}
}
\seealso{
\code{\link{modelExpandParam}, \link{scg}, \link{conjgrad}}.
}
\examples{
## missing
}
\keyword{model}
