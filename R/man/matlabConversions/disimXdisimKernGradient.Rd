\name{disimXdisimKernGradient}
\Rdversion{1.0}
\alias{disimXdisimKernGradient}
\title{Compute a cross gradient between two DISIM kernels.}
\description{
  Compute a cross gradient between two DISIM kernels.
}
\usage{
values <- disimXdisimKernGradient(disimKern1, disimKern2, t,
                                  covGrad)
values <- disimXdisimKernGradient(disimKern1, disimKern2, t1, t2,
                                  covGrad)
}
\arguments{
  \item{disimKern1}{the kernel structure associated with the first DISIM kernel.}
  \item{disimKern2}{the kernel structure associated with the second DISIM kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{covGrad}{gradient of the objective function with respect to the elements of the cross kernel matrix.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{g1, g2 <- disimXdisimKernGradient(disimKern1, disimKern2, t,
                                  covGrad)}
  computes cross gradient of parameters of a cross kernel between two disim kernels for the multiple output kernel.

  \code{g1, g2 <- disimXdisimKernGradient(disimKern1, disimKern2, t1, t2,
                                  covGrad)}
  computes cross kernel terms between two DISIM kernels for the multiple output kernel.

}
\value{
  \item{g1}{gradient of the parameters of the first kernel, for ordering see disimKernExtractParam.}
  \item{g2}{gradient of the parameters of the second kernel, for ordering see disimKernExtractParam.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{disimKernParamInit}, \link{disimKernExtractParam}}.
}
\examples{
## missing
}
\keyword{model}
