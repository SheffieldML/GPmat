\name{modelObjective}
\Rdversion{1.0}
\alias{modelObjective}
\title{Objective function to minimise for given model.}
\description{
  gives the objective function for a model. By default it is the negative log likelihood.
}
\usage{
err <- modelObjective(params, model, P1, P2, P3 ...)
}
\arguments{
  \item{params}{parameter vector to evaluate at.}
  \item{model}{model structure to optimise.}
  \item{P1, P2, P3 ...}{optional additional arguments.}
}
\value{
  \item{err}{the error function to be minimised. }
}
\seealso{
\code{\link{modelLogLikelihood}, \link{modelGradient}, \link{modelOptimise}}.
}
\examples{
## missing
}
\keyword{model}
