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
\alias{gammaPriorGradient}
\alias{invgammaPriorGradient}
\alias{kernPriorGradient}
\alias{priorGradient}
\alias{gammaPriorLogProb}
\alias{invgammaPriorLogProb}
\alias{kernPriorLogProb}
\alias{priorLogProb}
\title{Model log-likelihood/objective error function and its gradient.}
\description{
  \code{modeGradient} gives the gradient of the objective function for a
  model. By default the objective function (\code{modelObjective}) is
  a negative log likelihood (\code{modelLogLikelihood}).
}
\usage{
  modelObjective(params, model, ...)
  modelLogLikelihood(model)
  modelGradient(params, model, ...)
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
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # The probe identifier for TF 'twi'
  twi <- "143396_at"
  # The probe identifier for the target gene
  targetProbe <- "152715_at"

  # Create the model but do not optimise
  model <- GPLearn(drosophila_gpsim_fragment,
                   TF=twi, targets=targetProbe,
                   useGpdisim=TRUE, quiet=TRUE,
                   dontOptimise=TRUE)
  params <- modelExtractParam(model, only.values=FALSE)
  ll <- modelLogLikelihood(model)
  paramValues <- modelExtractParam(model, only.values=TRUE)
  modelGradient(paramValues, model)
}
\keyword{model}
