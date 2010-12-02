\name{gpBlockIndices}
\Rdversion{1.0}
\alias{gpBlockIndices}
\title{Return indices of given block.}
\description{
  returns the indices of a given block for the PITC approximation.
}
\usage{
indices <- gpBlockIndices(model, blockNum)
}
\arguments{
  \item{model}{the model for which the indices are being computed.}
  \item{blockNum}{the block number for which the indices are required.}
}
\value{
  \item{indices}{the data indices associated with given block.}
}
\seealso{
\code{\link{gpComputeAlpha}, \link{gpCovGrads}, \link{gpLogLikeGradients}, \link{gpLogLikelihood}, \link{gpUpdateAD}}.
}
\examples{
## missing
}
\keyword{model}
