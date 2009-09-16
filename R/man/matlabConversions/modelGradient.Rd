\name{modelGradient}
\Rdversion{1.0}
\alias{modelGradient}
\title{Gradient of error function to minimise for given model.}
\description{
  gives the gradient of the objective function for a model. By default the objective function is a negative log likelihood.
}
\usage{
g <- modelGradient(params, model, P1, P2, P3 ...)
}
\arguments{
  \item{params}{parameter vector to evaluate at.}
  \item{model}{model structure to optimise.}
  \item{P1, P2, P3 ...}{optional additional arguments.}
}
\value{
  \item{g}{the gradient of the error function to be minimised. }
}
\seealso{
\code{\link{modelLogLikeGradient}, \link{modelObjective}, \link{modelOptimise}}.
}
\examples{
## missing
}
\keyword{model}
