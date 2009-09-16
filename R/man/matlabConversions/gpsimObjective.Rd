\name{gpsimObjective}
\Rdversion{1.0}
\alias{gpsimObjective}
\title{Wrapper function for GPSIM objective.}
\description{
  returns the negative log likelihood of a Gaussian process model for single input motifs given the model structure and a vector parameters. This allows the use of NETLAB minimisation functions to find the model parameters.
}
\usage{
f <- gpsimObjective(params, model)
}
\arguments{
  \item{params}{the parameters of the model for which the objective will be evaluated.}
  \item{model}{the model structure for which the objective will be evaluated.}
}
\value{
  \item{f}{the negative log likelihood of the GPSIM model.}
}
\seealso{
\code{\link{scg}, \link{conjgrad}, \link{gpsimCreate}, \link{gpsimGradient}, \link{gpsimLogLikelihood}, \link{gpsimOptimise}}.
}
\examples{
## missing
}
\keyword{model}
