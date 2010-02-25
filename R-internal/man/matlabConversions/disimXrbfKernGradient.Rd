\name{disimXrbfKernGradient}
\Rdversion{1.0}
\alias{disimXrbfKernGradient}
\title{Compute gradient between the DISIM and RBF kernels.}
\description{
  Compute gradient between the DISIM and RBF kernels.
}
\usage{
values <- disimXrbfKernGradient(disimKern, rbfKern, t)
values <- disimXrbfKernGradient(disimKern, rbfKern, t1, t2)
}
\arguments{
  \item{disimKern}{the kernel structure associated with the DISIM kernel.}
  \item{rbfKern}{the kernel structure associated with the RBF kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{g1, g2 <- disimXrbfKernGradient(disimKern, rbfKern, t)}
  computes the gradient of an objective function with respect to cross kernel terms between DISIM and RBF kernels for the multiple output kernel.

  \code{g1, g2 <- disimXrbfKernGradient(disimKern, rbfKern, t1, t2)}
  computes the gradient of an objective function with respect to cross kernel terms between DISIM and RBF kernels for the multiple output kernel.

}
\value{
  \item{g1}{gradient of objective function with respect to kernel parameters of DISIM kernel.}
  \item{g2}{gradient of objective function with respect to kernel parameters of RBF kernel.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{disimKernParamInit}, \link{rbfKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
