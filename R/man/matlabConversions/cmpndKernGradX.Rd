\name{cmpndKernGradX}
\Rdversion{1.0}
\alias{cmpndKernGradX}
\title{Gradient of CMPND kernel with respect to a point x.}
\description{
  Gradient of CMPND kernel with respect to a point x.
}
\usage{
g <- cmpndKernGradX(kern, x)
g <- cmpndKernGradX(kern, x1, x2)
}
\arguments{
  \item{kern}{kernel structure for which gradients are being computed.}
  \item{x}{locations against which gradients are being computed.}
  \item{x1}{row locations against which gradients are being computed.}
  \item{x2}{column locations against which gradients are being computed.}
}
\details{
  \code{g <- cmpndKernGradX(kern, x)}
  computes the gradient of the compound kernel with respect to the input positions.

  \code{g <- cmpndKernGradX(kern, x1, x2)}
  computes the gradident of the compound kernel with respect to the input positions where both the row positions and column positions are provided separately.

}
\value{
  \item{g}{the returned gradients. The gradients are returned in a matrix which is numData x numInputs x numData. Where numData is the number of data points and numInputs is the number of input dimensions in X.}
}
\seealso{
\code{\link{}}.
}
\examples{
## missing
}
\keyword{model}
