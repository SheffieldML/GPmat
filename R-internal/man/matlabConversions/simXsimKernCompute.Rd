\name{simXsimKernCompute}
\Rdversion{1.0}
\alias{simXsimKernCompute}
\title{Compute a cross kernel between two SIM kernels.}
\description{
  Compute a cross kernel between two SIM kernels.
}
\usage{
values <- simXsimKernCompute(simKern1, simKern2, t)
values <- simXsimKernCompute(simKern1, simKern2, t1, t2)
}
\arguments{
  \item{simKern1}{the kernel structure associated with the first SIM kernel.}
  \item{simKern2}{the kernel structure associated with the second SIM kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{k, sK <- simXsimKernCompute(simKern1, simKern2, t)}
  computes cross kernel terms between two SIM kernels for the multiple output kernel.

  \code{k, sK <- simXsimKernCompute(simKern1, simKern2, t1, t2)}
  computes cross kernel terms between two SIM kernels for the multiple output kernel.

}
\value{
  \item{k}{block of values from kernel matrix.}
  \item{sK}{normalised matrix (i.e. unscaled version of K not multiplied by sqrt(simKern1.variance), sqrt(simKern2.variance) and sqrt(2/simKern.inverseWidth) in the case of the un-normalised kernel).}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{simKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
