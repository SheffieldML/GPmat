\name{whiteXwhiteKernGradient}
\Rdversion{1.0}
\alias{whiteXwhiteKernGradient}
\title{Compute a cross gradient between two WHITE kernels.}
\description{
  Compute a cross gradient between two WHITE kernels.
}
\usage{
values <- whiteXwhiteKernGradient(whiteKern1, whiteKern2, x,
                                  covGrad)
values <- whiteXwhiteKernGradient(whiteKern1, whiteKern2, x1, x2,
                                  covGrad)
}
\arguments{
  \item{whiteKern1}{the kernel structure associated with the first white noise kernel.}
  \item{whiteKern2}{the kernel structure associated with the second white noise kernel.}
  \item{x}{inputs for which kernel is to be computed.}
  \item{covGrad}{gradient of the objective function with respect to the elements of the cross kernel matrix.}
  \item{x1}{row inputs for which kernel is to be computed.}
  \item{x2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{g1, g2 <- whiteXwhiteKernGradient(whiteKern1, whiteKern2, x,
                                  covGrad)}
  computes cross gradient of parameters of a cross kernel between two white noise kernels for the multiple output kernel.

  \code{g1, g2 <- whiteXwhiteKernGradient(whiteKern1, whiteKern2, x1, x2,
                                  covGrad)}
  computes cross kernel terms between two white noise kernels for the multiple output kernel.

}
\value{
  \item{g1}{gradient of the parameters of the first kernel, for ordering see whiteKernExtractParam.}
  \item{g2}{gradient of the parameters of the second kernel, for ordering see whiteKernExtractParam.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{whiteKernParamInit}, \link{whiteKernExtractParam}}.
}
\examples{
## missing
}
\keyword{model}
