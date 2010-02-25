\name{disimXsimKernGradient}
\Rdversion{1.0}
\alias{disimXsimKernGradient}
\title{Compute gradient between the DISIM and SIM kernels.}
\description{
  Compute gradient between the DISIM and SIM kernels.
}
\usage{
values <- disimXsimKernGradient(disimKern, simKern, t)
values <- disimXsimKernGradient(disimKern, simKern, t1, t2)
}
\arguments{
  \item{disimKern}{the kernel structure associated with the DISIM kernel.}
  \item{simKern}{the kernel structure associated with the SIM kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{g1, g2 <- disimXsimKernGradient(disimKern, simKern, t)}
  computes the gradient of an objective function with respect to cross kernel terms between DISIM and SIM kernels for the multiple output kernel.

  \code{g1, g2 <- disimXsimKernGradient(disimKern, simKern, t1, t2)}
  computes the gradient of an objective function with respect to cross kernel terms between DISIM and SIM kernels for the multiple output kernel.

}
\value{
  \item{g1}{gradient of objective function with respect to kernel parameters of DISIM kernel.}
  \item{g2}{gradient of objective function with respect to kernel parameters of SIM kernel.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{simKernParamInit}, \link{rbfKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
