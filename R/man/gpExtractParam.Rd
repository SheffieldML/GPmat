\name{gpExtractParam}
\Rdversion{1.0}
\alias{gpExtractParam}
\title{Extract a parameter vector from a GP model.}
\description{
  does the same as above, but also returns parameter names.
}
\usage{
gpExtractParam(model, only.values=TRUE, ...)
}
\arguments{
  \item{model}{the model structure containing the information about the model.}
  \item{only.values}{(logical) do not return parameter names.}
  \item{...}{optional additional arguments.}
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
