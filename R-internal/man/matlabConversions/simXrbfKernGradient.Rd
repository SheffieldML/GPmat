\name{simXrbfKernGradient}
\Rdversion{1.0}
\alias{simXrbfKernGradient}
\title{Compute gradient between the SIM and RBF kernels.}
\description{
  Compute gradient between the SIM and RBF kernels.
}
\usage{
values <- simXrbfKernGradient(simKern, rbfKern, t)
values <- simXrbfKernGradient(simKern, rbfKern, t1, t2)
}
\arguments{
  \item{simKern}{the kernel structure associated with the SIM kernel.}
  \item{rbfKern}{the kernel structure associated with the RBF kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{g1, g2 <- simXrbfKernGradient(simKern, rbfKern, t)}
  computes the gradient of an objective function with respect to cross kernel terms between SIM and RBF kernels for the multiple output kernel.

  \code{g1, g2 <- simXrbfKernGradient(simKern, rbfKern, t1, t2)}
  computes the gradient of an objective function with respect to cross kernel terms between SIM and RBF kernels for the multiple output kernel.

}
\value{
  \item{g1}{gradient of objective function with respect to kernel parameters of SIM kernel.}
  \item{g2}{gradient of objective function with respect to kernel parameters of RBF kernel.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{simKernParamInit}, \link{rbfKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
