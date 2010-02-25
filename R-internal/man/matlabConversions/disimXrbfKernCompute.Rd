\name{disimXrbfKernCompute}
\Rdversion{1.0}
\alias{disimXrbfKernCompute}
\title{Compute a cross kernel between the DISIM and RBF kernels.}
\description{
  Compute a cross kernel between the DISIM and RBF kernels.
}
\usage{
k <- disimXrbfKernCompute(disimKern, rbfKern, t)
k <- disimXrbfKernCompute(disimKern, rbfKern, t1, t2)
}
\arguments{
  \item{disimKern}{the kernel structure associated with the DISIM kernel.}
  \item{rbfKern}{the kernel structure associated with the RBF kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{k <- disimXrbfKernCompute(disimKern, rbfKern, t)}
  computes cross kernel terms between DISIM and RBF kernels for the multiple output kernel.

  \code{k <- disimXrbfKernCompute(disimKern, rbfKern, t1, t2)}
  computes cross kernel terms between DISIM and RBF kernels for the multiple output kernel.

}
\value{
  \item{k}{block of values from kernel matrix.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{disimKernParamInit}, \link{rbfKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
