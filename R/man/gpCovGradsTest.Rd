\name{gpCovGradsTest}
\Rdversion{1.0}
\alias{gpCovGradsTest}
\title{Test the gradients of the likelihood wrt the covariance.}
\description{
  tests the gradients of the covariance to ensure they are correct.
}
\usage{
model <- gpCovGradsTest(model)
}
\arguments{
  \item{model}{the model to be tested.}
}
\value{
  \item{model}{the model that was tested.}
}
\seealso{
\code{\link{gpCreate}, \link{gpCovGrads}}.
}
\examples{
## missing
}
\keyword{model}
