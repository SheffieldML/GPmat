\name{gpTest}
\Rdversion{1.0}
\alias{gpTest}
\title{Test the gradients of the gpCovGrads function and the gp models.}
\description{
  runs some tests on the GP code to test that it is working.
}
\usage{
gpTest(q=2, d=3, N=10, k=5)
}
\arguments{
  \item{q}{input data dimension.}
  \item{d}{the number of processes (i.e. output data dimension).}
  \item{N}{the number of datapoints.}
  \item{k}{the number of inducing variables.}
}
\value{
  \item{model}{a cell array of models used for testing.}
}
% \seealso{
% \code{\link{modelTest}}.
% }
\examples{
## missing
}
\keyword{model}
