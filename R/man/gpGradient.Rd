\name{gpGradient}
\Rdversion{1.0}
\alias{gpGradient}
\title{Gradient wrapper for a GP model.}
\description{
  wraps the log likelihood gradient function to return the gradient of the negative of the log likelihood. This can then be used in, for example, NETLAB, minimisation tools.
}
\usage{
g <- gpGradient(params, model)
}
\arguments{
  \item{params}{the parameters of the model.}
  \item{model}{the model for which gradients will be computed.}
}
\value{
  \item{g}{the returned gradient of the negative log likelihood for the given parameters.}
}
\seealso{
  \code{\link{gpCreate}, \link{gpGradient}, \link{gpLogLikeGradients}, \link{gpOptimise}}.
}
\examples{
## missing
}
\keyword{model}
