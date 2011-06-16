\name{gpOut}
\Rdversion{1.0}
\alias{gpOut}
\title{Evaluate the output of an Gaussian process model.}
\description{
  evaluates the output of a given Gaussian process model.
}
\usage{
gpOut(model, x)
}
\arguments{
  \item{model}{the model for which the output is being evaluated.}
  \item{x}{the input position for which the output is required.}
}
\value{
  \item{y}{the output of the GP model. The function checks if there is a noise model associated with the GP, if there is, it is used, otherwise the mean of gpPosteriorMeanVar is returned.}
}
\seealso{
\code{\link{gpCreate}, \link{gpPosteriorMeanVar}}.
}
\examples{
## missing
}
\keyword{model}
