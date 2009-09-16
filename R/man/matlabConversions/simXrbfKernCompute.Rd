\name{simXrbfKernCompute}
\Rdversion{1.0}
\alias{simXrbfKernCompute}
\title{Compute a cross kernel between the SIM and RBF kernels.}
\description{
  Compute a cross kernel between the SIM and RBF kernels.
}
\usage{
values <- simXrbfKernCompute(simKern, rbfKern, t)
values <- simXrbfKernCompute(simKern, rbfKern, t1, t2)
}
\arguments{
  \item{simKern}{the kernel structure associated with the SIM kernel.}
  \item{rbfKern}{the kernel structure associated with the RBF kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{K, sK <- simXrbfKernCompute(simKern, rbfKern, t)}
  computes cross kernel terms between SIM and RBF kernels for the multiple output kernel.

  \code{K, sK <- simXrbfKernCompute(simKern, rbfKern, t1, t2)}
  computes cross kernel terms between SIM and RBF kernels for the multiple output kernel.

}
\value{
  \item{K}{block of values from kernel matrix.}
  \item{sK}{normalised matrix (i.e. unscaled version of K not multiplied by sqrt(kern.variance) and sqrt(2/simKern.inverseWidth) in the case of the un-normalised kernel).}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{simKernParamInit}, \link{rbfKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
