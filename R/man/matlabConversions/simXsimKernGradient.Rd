\name{simXsimKernGradient}
\Rdversion{1.0}
\alias{simXsimKernGradient}
\title{Compute a cross gradient between two SIM kernels.}
\description{
  Compute a cross gradient between two SIM kernels.
}
\usage{
values <- simXsimKernGradient(simKern1, simKern2, t, covGrad)
values <- simXsimKernGradient(simKern1, simKern2, t1, t2,
                              covGrad)
}
\arguments{
  \item{simKern1}{the kernel structure associated with the first SIM kernel.}
  \item{simKern2}{the kernel structure associated with the second SIM kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{covGrad}{gradient of the objective function with respect to the elements of the cross kernel matrix.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{g1, g2 <- simXsimKernGradient(simKern1, simKern2, t, covGrad)}
  computes cross gradient of parameters of a cross kernel between two sim kernels for the multiple output kernel.

  \code{g1, g2 <- simXsimKernGradient(simKern1, simKern2, t1, t2,
                              covGrad)}
  computes cross kernel terms between two SIM kernels for the multiple output kernel.

}
\value{
  \item{g1}{gradient of the parameters of the first kernel, for ordering see simKernExtractParam.}
  \item{g2}{gradient of the parameters of the second kernel, for ordering see simKernExtractParam.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{simKernParamInit}, \link{simKernExtractParam}}.
}
\examples{
## missing
}
\keyword{model}
