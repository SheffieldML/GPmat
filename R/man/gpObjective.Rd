\name{gpObjective}
\Rdversion{1.0}
\alias{gpObjective}
\title{Wrapper function for GP objective.}
\description{
  returns the negative log likelihood of a Gaussian process model given the model structure and a vector of parameters. This allows the use of NETLAB minimisation functions to find the model parameters.
}
\usage{
f <- gpObjective(params, model)
}
\arguments{
  \item{params}{the parameters of the model for which the objective will be evaluated.}
  \item{model}{the model structure for which the objective will be evaluated.}
}
\value{
  \item{f}{the negative log likelihood of the GP model.}
}
\seealso{
\code{
  \link{gpCreate}, \link{gpGradient}, \link{gpLogLikelihood}, \link{gpOptimise}}.
}
\examples{
## missing
}
\keyword{model}
