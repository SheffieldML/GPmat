\name{gpComputeAlpha}
\Rdversion{1.0}
\alias{gpComputeAlpha}
\title{Update the vector `alpha' for computing posterior mean quickly.}
\description{
  updates the vectors that are known as `alpha' in the support vector machine, in other words invK*y, where y is the target values.
}
\usage{
model <- gpComputeAlpha(model, m)
}
\arguments{
  \item{model}{the model for which the alphas are going to be updated.}
  \item{m}{the values of m for which the updates will be made.}
}
\value{
  \item{model}{the model with the updated alphas.}
}
\seealso{
\code{\link{gpCreate}, \link{gpUpdateAD}, \link{gpUpdateKernels}}.
}
\examples{
## missing
}
\keyword{model}
