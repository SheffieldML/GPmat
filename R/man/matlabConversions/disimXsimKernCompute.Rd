\name{disimXsimKernCompute}
\Rdversion{1.0}
\alias{disimXsimKernCompute}
\title{Compute a cross kernel between DISIM and SIM kernels.}
\description{
  Compute a cross kernel between DISIM and SIM kernels.
}
\usage{
k <- disimXsimKernCompute(disimKern, simKern, t)
k <- disimXsimKernCompute(disimKern, simKern, t1, t2)
}
\arguments{
  \item{disimKern}{the kernel structure associated with the DISIM kernel.}
  \item{simKern}{the kernel structure associated with the SIM kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{k <- disimXsimKernCompute(disimKern, simKern, t)}
  computes cross kernel terms between DISIM and SIM kernels for the multiple output kernel.

  \code{k <- disimXsimKernCompute(disimKern, simKern, t1, t2)}
  computes cross kernel terms between DISIM and SIM kernels for the multiple output kernel.

}
\value{
  \item{k}{block of values from kernel matrix.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{disimKernParamInit}, \link{simKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
