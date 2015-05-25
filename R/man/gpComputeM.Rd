\name{gpComputeM}
\Rdversion{1.0}
\alias{gpComputeM}
\title{Compute the matrix m given the model.}
\description{
  computes the matrix m (the scaled, bias and mean function removed matrix of the targets), given the model.
}
\usage{
  gpComputeM(model)
}
\arguments{
  \item{model}{the model for which the values are to be computed.}
}
\value{
  \item{m}{the scaled, bias and mean function removed values.}
}
\seealso{
\code{\link{gpCreate}, \link{gpComputeAlpha}, \link{gpUpdateAD}}.
}
\examples{
## missing
}
\keyword{model}
