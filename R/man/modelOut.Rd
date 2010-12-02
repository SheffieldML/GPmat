\name{modelOut}
\Rdversion{1.0}
\alias{modelOut}
\title{Give the output of a model for given X.}
\description{
  Give the output of a model for given X.
}
\usage{
Y <- modelOut(model, X)
values <- modelOut(model, X)
}
\arguments{
  \item{model}{structure specifying the model.}
  \item{X}{input location(s) for which output is to be computed.}
}
\details{
  \code{Y <- modelOut(model, X)}
  gives the output of the model for a given input X. For latent variable models it gives a position in data space given a position in latent space.

  \code{Phi, Y <- modelOut(model, X)}
  gives the output of the model for a given input X. For latent variable models it gives a position in data space given a position in latent space.

}
\value{
  \item{Y}{output location(s) corresponding to given input locations.}
  \item{Phi}{output basis function(s) corresponding to given input}
}
\seealso{
\code{\link{modelCreate}}.
}
\examples{
## missing
}
\keyword{model}
