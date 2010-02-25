\name{whiteXwhiteKernCompute}
\Rdversion{1.0}
\alias{whiteXwhiteKernCompute}
\title{Compute a cross kernel between two WHITE kernels.}
\description{
  Compute a cross kernel between two WHITE kernels.
}
\usage{
k <- whiteXwhiteKernCompute(whiteKern1, whiteKern2, x)
k <- whiteXwhiteKernCompute(whiteKern1, whiteKern2, x1, x2)
}
\arguments{
  \item{whiteKern1}{the kernel structure associated with the first white noise kernel.}
  \item{whiteKern2}{the kernel structure associated with the second white noise kernel.}
  \item{x}{inputs for which kernel is to be computed.}
  \item{x1}{row inputs for which kernel is to be computed.}
  \item{x2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{k <- whiteXwhiteKernCompute(whiteKern1, whiteKern2, x)}
  computes cross kernel terms between two white noise kernels for the multiple output kernel.

  \code{k <- whiteXwhiteKernCompute(whiteKern1, whiteKern2, x1, x2)}
  computes cross kernel terms between two WHITE kernels for the multiple output kernel.

}
\value{
  \item{k}{block of values from kernel matrix.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{whiteKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
