\name{gpPosteriorMeanVar}
\Rdversion{1.0}
\alias{gpPosteriorMeanVar}
\title{Mean and variances of the posterior at points given by X.}
\description{
  returns the posterior mean and variance for a given set of points.
}
\usage{
values <- gpPosteriorMeanVar(model, x)
}
\arguments{
  \item{model}{the model for which the posterior will be computed.}
  \item{x}{the input positions for which the posterior will be computed.}
}
\value{
  \item{mu}{the mean of the posterior distribution.}
  \item{sigma}{the variances of the posterior distributions.}
}
\seealso{
\code{
  \link{gpCreate}}.
}
\examples{
## missing
}
\keyword{model}
