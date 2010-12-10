\name{gpExtractParam}
\Rdversion{1.0}
\alias{gpExtractParam}
\title{Extract a parameter vector from a GP model.}
\description{
  does the same as above, but also returns parameter names.
}
\usage{
values <- gpExtractParam(model, model)
}
\arguments{
  \item{model}{the model structure containing the information about the model.}
}
\value{
  \item{params}{a vector of parameters from the model.}
  \item{names}{cell array of parameter names.}
}
\seealso{
\code{\link{gpCreate}, \link{gpExpandParam}, \link{modelExtractParam}}.
}
\examples{
## missing
}
\keyword{model}
