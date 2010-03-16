\name{modelGradient}
\Rdversion{1.0}
\alias{modelGradient}
\alias{modelObjective}
\alias{modelLogLikelihood}
\alias{cgpdisimGradient}
\alias{cgpdisimLogLikeGradients}
\alias{cgpdisimObjective}
\alias{cgpdisimLogLikelihood}
\alias{cgpsimGradient}
\alias{cgpsimLogLikeGradients}
\alias{cgpsimObjective}
\alias{cgpsimLogLikelihood}
\alias{gpdisimGradient}
\alias{gpdisimLogLikeGradients}
\alias{gpdisimObjective}
\alias{gpdisimLogLikelihood}
\alias{gpsimGradient}
\alias{gpsimLogLikeGradients}
\alias{gpsimObjective}
\alias{gpsimLogLikelihood}
\title{Model log-likelihood/objective error function and its gradient.}
\description{
  \code{modeGradient} gives the gradient of the objective function for a
  model. By default the objective function (\code{modelObjective}) is
  a negative log likelihood (\code{modelLogLikelihood}).
}
\usage{
v <- modelObjective(model)
ll <- modelLogLikelihood(model)
g <- modelGradient(params, model, ...)
}
\arguments{
  \item{params}{parameter vector to evaluate at.}
  \item{model}{model structure.}
  \item{...}{optional additional arguments.}
}
\value{
  \item{g}{the gradient of the error function to be minimised. }
  \item{v}{the objective function value (lower is better).}
  \item{ll}{the log-likelihood value.}
}
\seealso{
\code{\link{modelOptimise}}.
}
\examples{
## missing
}
\keyword{model}
